import os,time,sys, pwd, random, datetime, shutil
todaysdate = datetime.date.today()

uname = pwd.getpwuid(os.getuid()).pw_name

def makeSubmitScript(toyDir, seed, n):
	initDir=os.getcwd();
	dirname=toyDir
	os.chdir(dirname)
	#os.system("echo '#!/bin/bash' >> BatchSubmit.sh");
#	os.system("echo '#PBS -l pmem=3000mb' >> BatchSubmit.sh");
	#os.system("echo '#PBS -l ncpus=1' >> BatchSubmit.sh");
	#os.system("echo '. /data/lhcb/sw/scripts/lbsetup-cvmfs-osagnostic.sh\n' >> BatchSubmit.sh");
	#os.system("echo '. LbLogin.sh -c x86_64-slc6-gcc49-opt\n' >> BatchSubmit.sh");
	#os.system("echo '. SetupProject.sh Gaudi v26r1 ROOT \n ' >> BatchSubmit.sh")

    os.system("echo '#!/bin/bash' >> BatchSubmit.sh");
    os.system("echo '#PBS -l cput=6:00:00' >> BatchSubmit.sh");
    os.system("echo '#PBS -l walltime=6:00:00' >> BatchSubmit.sh");
    os.system("echo 'source /data/lhcb/sw/scripts/lbsetup-cvmfs.sh' >> BatchSubmit.sh");
    os.system("echo '. SetupProject.sh Gaudi ROOT' >> BatchSubmit.sh");

	os.system("echo 'cd "+dirname+"/../Bu2D0Kstar_hh_MassFit_copy/' >> BatchSubmit.sh") ## FIX
	os.system("echo './bin/run ../toySeed_"+seed+"/GeneralSettings.txt > ../toySeed_"+seed+"/runOutput.txt'>> BatchSubmit.sh");
    os.system("echo 'bin/run -m Ds2PhiPi:Ds2KKPi:Ds2PiPiPi:Ds2KPiPi -d s21:s21r1:s24 -S -M -c DsD0:DsPhi:DsPhiSide -Z -t "+n+" -T toysDir/ -s "+seed+"'>> BatchSubmit.sh");
	os.chdir(initDir)


basedir = "/data/lhcb/users/" + uname + "/B2DKstar/ToyStudies"
if not os.path.exists(basedir):
    print 'Creating', basedir
    os.mkdir(basedir)
    os.chmod(basedir, 0750)

desc = raw_input("Enter a description of the toys you are running (today's date will automatically be included): ")

if desc == '':
    print 'You entered a blank string so I am using "toys"'
    desc = 'toys'

toyDir = basedir + '/' + str(todaysdate) + "_" + desc
if os.path.exists(toyDir):
    print
    print '*** WARNING: directory', toyDir, 'already exists! ***'
    decision = raw_input('*** Input Y to overwrite, anything else to abort: ').lower()
    if decision != 'y':
        sys.exit()
    else:
        shutil.rmtree(toyDir)
        print 'Directory has been removed and will be recreated'
print 'Creating', toyDir
os.mkdir(toyDir)
print
# Get seed
seedInput = int(raw_input("Give random seed: "))

numBatchJobs = int(raw_input("Num batch jobs: "))
numMCStudyToysPerBatchJob = int(raw_input("Num MCStudyToys per batch job (>1): "))

##	Tar up and backup the source code there:
##os.system('tar -czf gammaFitToys/batchToyScriptsAndOutput/'+toyDir+"/codeSnapshot.tgz Bu2D0H_KSHH_BinnedFit2013/");
print 'Copying the code into the new toy directory.' 
# Make a copy of the code there ready to run - then don't have to wait for each toy to complete
os.system('mkdir '+toyDir+'/Bu2D0Kstar_hh_MassFit_copy')
os.system("cp -r Bu2D0Kstar_hh_MassFit/* "+toyDir+"/Bu2D0Kstar_hh_MassFit_copy")
	# remove unnecessary files from the copy directory%
print "I'm about to clear out the copied code directory ("+toyDir+"/Bu2D0Kstar_hh_MassFit_copy/). Please confirm the following 'rm' commands if they seem sensible:"
print "    rm -f "+toyDir+"/Bu2D0Kstar_hh_MassFit_copy/*.o"
print "    rm -f "+toyDir+"/Bu2D0Kstar_hh_MassFit_copy/.svn"
print "    rm -f "+toyDir+"/Bu2D0Kstar_hh_MassFit_copy/figs/*"
print "    rm -f "+toyDir+"/Bu2D0Kstar_hh_MassFit_copy/figs/fits/*"
print "    rm -f "+toyDir+"/Bu2D0Kstar_hh_MassFit_copy/figs/residuals/*"

if(raw_input("Are you sure you want to continue? (y/n) ")!="y"):
	sys.exit();
os.system("rm -f "+toyDir+"/Bu2D0Kstar_hh_MassFit_copy/*.o")
os.system("rm -f "+toyDir+"/Bu2D0Kstar_hh_MassFit_copy/.svn")
os.system("rm -f "+toyDir+"/Bu2D0Kstar_hh_MassFit_copy/figs/*")
os.system("rm -f "+toyDir+"/Bu2D0Kstar_hh_MassFit_copy/figs/fits/*")
os.system("rm -f "+toyDir+"/Bu2D0Kstar_hh_MassFit_copy/figs/residuals/*")

# Start toys
jobCounter=1
for toySeed in range(seedInput,numBatchJobs*numMCStudyToysPerBatchJob+seedInput,numMCStudyToysPerBatchJob):
	# Make toy directory
	toyDirName=toyDir+"/toySeed_"+str(toySeed)+"/"
	os.system('mkdir '+toyDirName)

	# Create copy of fit code parameter file there and set the seed
        os.system('cp Bu2D0Kstar_hh_MassFit/Settings/GeneralSettings.txt '+toyDirName+'GeneralSettings.txt')
	f_fitParFile=open(toyDirName+'GeneralSettings.txt')
	f_newFitParFile=open(toyDirName+'GeneralSettings_new.txt','w')
	for line in f_fitParFile.readlines():
		if(line.find('doFit')==0):
			f_newFitParFile.write('doFit true\n')
		elif(line.find('drawProjections')==0):
			f_newFitParFile.write('drawProjections false\n')
		elif(line.find('startSeed')==0):
			f_newFitParFile.write('startSeed '+str(toySeed)+"\n")
		elif(line.find('nToys')==0):
			f_newFitParFile.write('nToys '+str(numMCStudyToysPerBatchJob)+'\n')
		elif(line.find('genToys')==0):
			f_newFitParFile.write('genToys true\n')
		elif(line.find('readToys')==0):
			f_newFitParFile.write('readToys false\n')
		elif(line.find('toyLocation')==0):
			f_newFitParFile.write('toyLocation ../toySeed_'+str(toySeed)+'/\n')
		elif(line.find('UNBLIND')==0):
			f_newFitParFile.write('UNBLIND true\n')
		else:	f_newFitParFile.write(line)
	f_fitParFile.close();
	f_newFitParFile.close();
	os.system('cd '+toyDirName+'; mv GeneralSettings_new.txt GeneralSettings.txt')

	# Make batch submit script0
	print toyDirName
	makeSubmitScript(toyDirName, str(toySeed))

# Submit script from that directory
	initDir=os.getcwd();
	os.chdir(toyDirName)
	print 'Submitting job '+str(jobCounter)+" with seed "+str(toySeed)
	time.sleep(0.1)
	os.system('qsub -lcput=1:59:59 -N Toy_'+str(toySeed)+' BatchSubmit.sh')
	#os.system('qsub -q testing -lcput=1:59:59 -N Toy_'+str(toySeed)+' BatchSubmit.sh')
#	os.system('. BatchSubmit.sh &')
	os.chdir(initDir)
	jobCounter+=1
