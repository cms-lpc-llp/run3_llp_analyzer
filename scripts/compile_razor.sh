#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
script_name="$(basename "${BASH_SOURCE[0]}")"
log_prefix="${script_name%.sh}"
container_launcher="/cvmfs/cms.cern.ch/common/cmssw-el9"

if [ -z "${RUN3RUNNER_IN_CONTAINER:-}" ]; then
  if [ "${RUN3_FORCE_NATIVE:-0}" = "1" ]; then
    echo "[INFO] RUN3_FORCE_NATIVE=1; running natively on host."
  else
    if [ ! -x "${container_launcher}" ]; then
      echo "[ERROR] Missing container launcher: ${container_launcher}" >&2
      echo "[ERROR] Set RUN3_FORCE_NATIVE=1 to skip container relaunch." >&2
      exit 1
    fi
    echo "[INFO] Relaunching inside el9 CMSSW container..."
    exec "${container_launcher}" --command-to-run \
      "bash -lc \"export RUN3RUNNER_IN_CONTAINER=1; cd \\\"${script_dir}\\\"; bash \\\"${script_name}\\\"\""
  fi
fi

analyzer_dir_default="$(cd "${script_dir}/.." && pwd)"
root_setup_default="/cvmfs/cms.cern.ch/el9_amd64_gcc12/lcg/root/6.30.07-024df6516c17fd2edef848a927a788f1/bin/thisroot.sh"

if [ -n "${RUN3_ANALYZER_DIR:-}" ]; then
  analyzer_dir="${RUN3_ANALYZER_DIR}"
elif [ -n "${CMSSW_BASE:-}" ]; then
  analyzer_dir="${CMSSW_BASE}/src/run3_llp_analyzer"
else
  analyzer_dir="${analyzer_dir_default}"
fi

if [ -n "${RUN3_CMSSW_SRC:-}" ]; then
  cmssw_src="${RUN3_CMSSW_SRC}"
elif [ -n "${CMSSW_BASE:-}" ]; then
  cmssw_src="${CMSSW_BASE}/src"
else
  cmssw_src="$(cd "${analyzer_dir}/.." && pwd)"
fi

if [ -n "${RUN3_ROOT_SETUP:-}" ]; then
  root_setup="${RUN3_ROOT_SETUP}"
else
  root_setup="${root_setup_default}"
fi

if [ -n "${RUN3_LOG_DIR:-}" ]; then
  log_dir="${RUN3_LOG_DIR}"
elif [ -n "${CMSSW_BASE:-}" ]; then
  log_dir="${CMSSW_BASE}/logs"
else
  log_dir="$(cd "${cmssw_src}/.." && pwd)/logs"
fi

mkdir -p "${log_dir}"
rm -f "${log_dir}/${log_prefix}.out" \
      "${log_dir}/${log_prefix}.err" \
      "${log_dir}/${log_prefix}.wrn" || true

out_log="${log_dir}/${log_prefix}.out"
wrn_log="${log_dir}/${log_prefix}.wrn"
err_log="${log_dir}/${log_prefix}.err"
stderr_tmp="$(mktemp "${log_dir}/${log_prefix}.stderr.XXXXXX")"

: >"${out_log}"
: >"${wrn_log}"
: >"${err_log}"

split_stderr_logs() {
  if [ -s "${stderr_tmp}" ]; then
    grep -i 'warning' "${stderr_tmp}" > "${wrn_log}" || true
    grep -iv 'warning' "${stderr_tmp}" > "${err_log}" || true
  else
    : >"${wrn_log}"
    : >"${err_log}"
  fi
}

finalize_logs() {
  local status=$?
  split_stderr_logs
  rm -f "${stderr_tmp}"
  return "${status}"
}

trap 'finalize_logs' EXIT

exec 3>&1
exec 1>>"${out_log}" \
  2> >(tee -a "${stderr_tmp}" >> "${out_log}")

step() {
  local stamp line
  stamp="$(date +'%F %T')"
  line="[$stamp] $1"
  printf '%s\n' "${line}" >> "${out_log}"
  if [ -w /dev/tty ]; then
    printf '%s\n' "${line}" > /dev/tty
  else
    printf '%s\n' "${line}" >&3
  fi
}

step "Using existing analyzer checkout (no clone/cmsrel/git prep)"
step "Analyzer directory: ${analyzer_dir}"
step "CMSSW src directory: ${cmssw_src}"

if [ ! -d "${analyzer_dir}" ]; then
  echo "[ERROR] Analyzer directory not found: ${analyzer_dir}" >&2
  echo "[ERROR] Set RUN3_ANALYZER_DIR or CMSSW_BASE to an existing checkout." >&2
  exit 1
fi

if [ ! -f "${analyzer_dir}/Makefile" ]; then
  echo "[ERROR] Missing Makefile in analyzer directory: ${analyzer_dir}" >&2
  exit 1
fi

if [ ! -d "${cmssw_src}" ]; then
  echo "[ERROR] CMSSW src directory not found: ${cmssw_src}" >&2
  echo "[ERROR] Set RUN3_CMSSW_SRC or CMSSW_BASE to an existing CMSSW area." >&2
  exit 1
fi

step "Activating CMSSW runtime and ROOT"
set +u
source /cvmfs/cms.cern.ch/cmsset_default.sh
set -u
export SCRAM_ARCH="${RUN3_SCRAM_ARCH:-el9_amd64_gcc12}"

cd "${cmssw_src}"
if ! eval "$(scramv1 runtime -sh)"; then
  echo "[ERROR] scram runtime activation failed in ${cmssw_src}" >&2
  exit 1
fi

if [ -f "${root_setup}" ]; then
  # Keep ROOT environment consistent with the EL9 build toolchain.
  source "${root_setup}"
else
  echo "[ERROR] ROOT setup script not found: ${root_setup}" >&2
  echo "[ERROR] Override with RUN3_ROOT_SETUP if your ROOT path differs." >&2
  exit 1
fi

cd "${analyzer_dir}"

step "Running make clean"
make clean

# Parallel build (dependency ordering is handled in the analyzer Makefile).
if command -v nproc >/dev/null 2>&1; then
  default_jobs="$(nproc)"
else
  default_jobs=8
fi
jobs="${RUN3_MAKE_JOBS:-${default_jobs}}"

case "${jobs}" in
  ''|*[!0-9]*)
    echo "[ERROR] RUN3_MAKE_JOBS must be a positive integer (got '${jobs}')." >&2
    exit 1
    ;;
  0)
    echo "[ERROR] RUN3_MAKE_JOBS must be greater than zero." >&2
    exit 1
    ;;
esac

step "Compiling/linking with make -s -j${jobs} all"
step "Progress tracked in ${log_prefix}.out"
if make -s -j"${jobs}" all; then
  step "Build succeeded"
  exit 0
else
  step "Build failed; see ${log_prefix}.err and ${log_prefix}.wrn"
  exit 1
fi
