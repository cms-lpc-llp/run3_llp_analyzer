#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  submit_ucsd.sh --input-list <path> [options]

Required:
  --input-list <path>            Path to input .txt list

Options:
  --cmssw-base <path>            CMSSW base (default: $CMSSW_BASE)
  --compile-script <relpath>     Relative to run3_llp_analyzer (default: scripts/compile_razor.sh)
  --skip-compile                 Skip compile step
  --analyzer <name>              Analyzer name (default: llp_MuonSystem_CA_mdsnano)
  --files-per-job <n>            Files per condor job (default: 1)
  --is-data <no|yes|--isData>    Data flag passed to runner (default: no)
  --analyzer-tag <tag>           Analyzer tag (default: Summer24)
  --output-directory <path>      Output destination (default: /ceph/cms/store/group/mds-ml/run3_analyzer_output/)
  --job-flavour <name>           HTCondor JobFlavour (default: espresso)
  --desired-site <name>          HTCondor site preference (default: T2_US_UCSD)
  --request-memory <mb>          RequestMemory (default: 1000)
  --request-cpus <n>             RequestCpus (default: 1)
  --request-disk <kb>            RequestDisk KB (default: 2000000)
  --singularity-image <img>      Singularity image (default: /cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel9)
  --x509-proxy <path>            Proxy path (default: $X509_USER_PROXY or /tmp/x509up_u81902)
  --test-mode                    Condor test mode: first non-empty input line only, 1 job
  --batch-name <name>            Condor batch name override
  --dry-run                      Write JDL and print submit command only
  -h, --help                     Show this help
EOF
}

run_cmd() {
  echo "+ $*"
  "$@"
}

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
analyzer_dir="$(cd "${script_dir}/.." && pwd)"
script_name="$(basename "${BASH_SOURCE[0]}")"
script_stem="${script_name%.*}"

input_list=""
cmssw_base="${CMSSW_BASE:-}"
compile_script="scripts/compile_razor.sh"
skip_compile=0
analyzer="llp_MuonSystem_CA_mdsnano"
files_per_job=1
is_data="no"
analyzer_tag="Summer24"
output_directory="/ceph/cms/store/group/mds-ml/run3_analyzer_output/"
job_flavour="espresso"
desired_site="T2_US_UCSD"
request_memory=1000
request_cpus=1
request_disk=2000000
singularity_image="/cvmfs/singularity.opensciencegrid.org/cmssw/cms:rhel9"
x509_proxy="${X509_USER_PROXY:-/tmp/x509up_u81902}"
test_mode=0
batch_name=""
dry_run=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-list) input_list="${2:-}"; shift 2 ;;
    --cmssw-base) cmssw_base="${2:-}"; shift 2 ;;
    --compile-script) compile_script="${2:-}"; shift 2 ;;
    --skip-compile) skip_compile=1; shift ;;
    --analyzer) analyzer="${2:-}"; shift 2 ;;
    --files-per-job) files_per_job="${2:-}"; shift 2 ;;
    --is-data) is_data="${2:-}"; shift 2 ;;
    --analyzer-tag) analyzer_tag="${2:-}"; shift 2 ;;
    --output-directory) output_directory="${2:-}"; shift 2 ;;
    --job-flavour) job_flavour="${2:-}"; shift 2 ;;
    --desired-site) desired_site="${2:-}"; shift 2 ;;
    --request-memory) request_memory="${2:-}"; shift 2 ;;
    --request-cpus) request_cpus="${2:-}"; shift 2 ;;
    --request-disk) request_disk="${2:-}"; shift 2 ;;
    --singularity-image) singularity_image="${2:-}"; shift 2 ;;
    --x509-proxy) x509_proxy="${2:-}"; shift 2 ;;
    --test-mode) test_mode=1; shift ;;
    --batch-name) batch_name="${2:-}"; shift 2 ;;
    --dry-run) dry_run=1; shift ;;
    -h|--help) usage; exit 0 ;;
    *)
      echo "Unknown option: $1" >&2
      usage
      exit 2
      ;;
  esac
done

if [[ -z "${input_list}" ]]; then
  echo "--input-list is required" >&2
  usage
  exit 2
fi
if [[ -z "${cmssw_base}" ]]; then
  echo "--cmssw-base not provided and CMSSW_BASE is unset" >&2
  exit 2
fi

input_list="$(readlink -f "${input_list}")"
cmssw_base="$(readlink -f "${cmssw_base}")"

if [[ ! -f "${input_list}" ]]; then
  echo "input list does not exist: ${input_list}" >&2
  exit 2
fi
if [[ ! -d "${cmssw_base}" ]]; then
  echo "CMSSW base does not exist: ${cmssw_base}" >&2
  exit 2
fi

submitter_log_dir="${cmssw_base}/logs"
condor_tmp_dir="${cmssw_base}/tmp/${script_stem}"
mkdir -p "${submitter_log_dir}" "${condor_tmp_dir}"
submitter_out_log="${submitter_log_dir}/${script_stem}.out"
submitter_err_log="${submitter_log_dir}/${script_stem}.err"
: > "${submitter_out_log}"
: > "${submitter_err_log}"
exec > >(tee -a "${submitter_out_log}") 2> >(tee -a "${submitter_err_log}" >&2)

echo "submitter out log: ${submitter_out_log}"
echo "submitter err log: ${submitter_err_log}"

if (( test_mode )); then
  test_list_dir="${condor_tmp_dir}/test_lists"
  mkdir -p "${test_list_dir}"
  first_line="$(awk 'NF {print; exit}' "${input_list}")"
  if [[ -z "${first_line}" ]]; then
    echo "input list has no non-empty lines: ${input_list}" >&2
    exit 2
  fi
  test_list="${test_list_dir}/$(basename "${input_list}" .txt)_condor_test.txt"
  printf "%s\n" "${first_line}" > "${test_list}"
  echo "test_mode input file: ${first_line}"
  echo "test_mode list path: ${test_list}"
  input_list="${test_list}"
  files_per_job=1
fi

sample="$(basename "${input_list}" .txt)"
nfiles="$(awk 'NF {c++} END {print c+0}' "${input_list}")"
if [[ "${nfiles}" -eq 0 ]]; then
  echo "input list has no non-empty lines: ${input_list}" >&2
  exit 2
fi
maxjob="$(( (nfiles + files_per_job - 1) / files_per_job ))"

echo "sample=${sample}"
echo "input_list=${input_list}"
echo "nfiles=${nfiles} filesPerJob=${files_per_job} maxjob=${maxjob}"
echo "cmssw_base=${cmssw_base}"
echo "output_directory=${output_directory}"
if (( test_mode )); then
  echo "mode=condor_test (single file, single job)"
fi

if (( ! skip_compile )); then
  compile_script_path="${analyzer_dir}/${compile_script}"
  if [[ ! -f "${compile_script_path}" ]]; then
    echo "compile script not found: ${compile_script_path}" >&2
    exit 2
  fi
  export CMSSW_BASE="${cmssw_base}"
  run_cmd bash "${compile_script_path}"
fi

timestamp="$(date +%Y%m%d_%H%M%S)"
tarball_dir="${condor_tmp_dir}/tarballs"
mkdir -p "${tarball_dir}"
cmssw_name="$(basename "${cmssw_base}")"
tarball_path="${tarball_dir}/${cmssw_name}_compiled_${timestamp}.tar.gz"
run_cmd tar \
  --exclude="${cmssw_name}/logs" \
  --exclude="${cmssw_name}/tmp/${script_stem}" \
  --exclude="${cmssw_name}/src/run3_llp_analyzer/scripts_condor/log" \
  --exclude="${cmssw_name}/src/run3_llp_analyzer/scripts_condor/tarballs" \
  -C "$(dirname "${cmssw_base}")" \
  -czf "${tarball_path}" \
  "${cmssw_name}"
echo "tarball=${tarball_path}"

submit_dir="${condor_tmp_dir}/submit"
condor_log_dir="${submitter_log_dir}/condor"
mkdir -p "${submit_dir}" "${condor_log_dir}"

jdl_file="${submit_dir}/ucsd_${sample}_${maxjob}.jdl"
executable="untar.sh"
run_script="${script_dir}/runAnalyzer_ucsd.sh"
run_wrapper_script="${script_dir}/untar.sh"

home_dir="${HOME:-}"
if [[ -n "${home_dir}" && "${home_dir}" != */ ]]; then
  home_dir="${home_dir}/"
fi
cmssw_base_hint="from_tarball:${cmssw_name}"

arguments="${analyzer} $(basename "${input_list}") ${is_data} ${files_per_job} \$(ProcId) ${maxjob} ${output_directory} ${analyzer_tag} ${cmssw_base_hint} ${home_dir} $(basename "${tarball_path}")"
transfer_input_files="${run_wrapper_script},${run_script},${input_list},${tarball_path}"

cat > "${jdl_file}" <<EOF
Universe = vanilla
Executable = ${executable}
Arguments = ${arguments}
Log = ${condor_log_dir}/ucsd_${sample}_Job\$(ProcId)_Of_${maxjob}_\$(Cluster).\$(Process).log
Output = ${condor_log_dir}/ucsd_${sample}_Job\$(ProcId)_Of_${maxjob}_\$(Cluster).\$(Process).out
Error = ${condor_log_dir}/ucsd_${sample}_Job\$(ProcId)_Of_${maxjob}_\$(Cluster).\$(Process).err
+JobFlavour = "${job_flavour}"
+DESIRED_Sites = "${desired_site}"
RequestMemory = ${request_memory}
RequestCpus = ${request_cpus}
RequestDisk = ${request_disk}
+RunAsOwner = True
+SingularityImage = "${singularity_image}"
+SingularityBindCVMFS = True
use_x509userproxy = True
x509userproxy = ${x509_proxy}
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
StreamOut = True
StreamErr = True
Environment = "RUN_DIAG=1"
Transfer_Input_Files = ${transfer_input_files}
Queue ${maxjob}
EOF

echo "wrote JDL: ${jdl_file}"

if [[ -z "${batch_name}" ]]; then
  if (( test_mode )); then
    batch_name="ucsd_test_${sample}"
  else
    batch_name="ucsd_${sample}"
  fi
fi

submit_cmd=(condor_submit "${jdl_file}" --batch-name "${batch_name}")
if (( dry_run )); then
  echo "+ ${submit_cmd[*]}  # dry-run (not executed)"
else
  run_cmd "${submit_cmd[@]}"
fi
