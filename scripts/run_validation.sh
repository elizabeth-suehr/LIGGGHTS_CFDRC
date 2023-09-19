./clear.sh
cd rod2
ln -s ../../src/lmp_auto
cd ..
cd rod6
ln -s ../../src/lmp_auto
cd ..
cd multisphere_13
ln -s ../../src/lmp_auto
cd ..
python3 validation.py