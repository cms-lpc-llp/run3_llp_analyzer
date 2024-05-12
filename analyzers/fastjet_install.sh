curr_dir=$(pwd);

echo "Downloading FastJet";
wget wget https://fastjet.fr/repo/fastjet-3.4.2.tar.gz;
tar -xzf fastjet-3.4.2.tar.gz;

echo "Installing FastJet";
cd fastjet-3.4.2;
./configure --prefix=$curr_dir/fastjet-install;
make;
make check;
make install;

cd ..;
rm -rf fastjet-3.4.2.tar.gz;
echo "Installed FastJet successfully";
echo "FastJet installed at $curr_dir/fastjet-install";

cd $curr_dir;