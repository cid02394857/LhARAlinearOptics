# UserDirectory

The idea of this UserDirectory is that you can copy it to a scratch area and use the code in the main repository
in a way that keeps your work separate from the main logic of the software.

## To set up and run:

The first thing that you should do is to follow the setup instructions in the main repository.

* You should already have created the virtual environment (venv).
* You need to activate the virtual environment (source venv/bin/activate).
* You should have already installed the files in "requirements.txt". If not do so ("pip install -r requirements.txt").
* If you are using Anaconda and have (base) activated I recommend disabling with "conda deactivate".
* If you are on a Mac and "echo $SHELL" displays "/bin/zsh" you should start a bash shell by running "bash".
* For workflows based on this UserDirectory it's not necessary to source the startup.bash script in the main repository.

Once you have done all that:

* Copy this directory to a scratch area, for example /home/user/scratch: "cp -r 31-UserDirectory /home/user/scratch".
* Change to that directory: "cd /home/user/scratch".
* Assuming that you had checked out the main repo to /home/user/Git/LhARAlinearOptics, source the startup.bash file as follows:
  "source user_startup.bash -p /home/user/Git/LhARAlinearOptics"

You should now hopefully be ready to start using the scripts in the UserDirectory.
As an example:

$LhARAOpticsPATH/03-Scripts/runBeamSim.py -b 11-Parameters/LIONBeamLine-Params-Gauss.csv -o 99-Scratch/test.dat -n 10000

(Here we use the $LhARAOpticsPATH environment variable to run a script in the main repository passing some sample parameters from the UserDirectory)

Next we can run:

03-Scripts/visualiseBeam.py -i 99-Scratch/test.dat -o 99-Scratch/test.pdf -n 1000

(Here we are running the "visualiseBeam.py" script in the UserDirectory itself)
