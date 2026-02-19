curr_dir=$(pwd);

echo "Downloading FastJet";
wget https://fastjet.fr/repo/fastjet-3.4.2.tar.gz;
tar -xzf fastjet-3.4.2.tar.gz;

echo "Installing FastJet";
cd fastjet-3.4.2;
if command -v nproc >/dev/null 2>&1; then
  default_jobs="$(nproc)"
else
  default_jobs=8
fi
jobs="${FASTJET_MAKE_JOBS:-${default_jobs}}"
./configure --prefix=$curr_dir/fastjet-install;
make -j"${jobs}";
make -j"${jobs}" check;
make install;

cd ..;
rm -rf fastjet-3.4.2.tar.gz;
echo "Installed FastJet successfully";
echo "FastJet installed at $curr_dir/fastjet-install";

cd $curr_dir;
