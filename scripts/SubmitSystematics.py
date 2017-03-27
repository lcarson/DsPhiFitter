import os,time,sys, pwd, random, datetime, shutil
todaysdate = datetime.date.today()

def makeSubmitScript(sysDir, seed, variableList, wall, brstring, inputMode,others):
    initDir=os.getcwd();
    dirname=sysDir
    os.chdir(dirname)
    os.system("echo '#!/bin/bash' >> BatchSubmit_"+seed+".sh");
    os.system("echo '#PBS -l cput="+wall+":00:00' >> BatchSubmit_"+seed+".sh");
    os.system("echo '#PBS -l walltime="+wall+":00:00' >> BatchSubmit_"+seed+".sh");
    os.system("echo 'source /data/lhcb/sw/scripts/lbsetup-cvmfs.sh' >> BatchSubmit_"+seed+".sh");
    os.system("echo '. SetupProject.sh Gaudi ROOT' >> BatchSubmit_"+seed+".sh");

    os.system("echo 'cd "+dirname+"/../' >> BatchSubmit_"+seed+".sh")
    if inputMode == 'DsPhi':
        os.system("echo 'bin/run -m Ds2PhiPi:Ds2KKPi:Ds2PiPiPi:Ds2KPiPi -L data/FullSel -d s21:s21r1:s24:s26 "+others+" -c DsD0:DsPhi:DsPhiSide "+brstring+" -v "+variableList+" -V systematicsDir -s "+seed+"'>> BatchSubmit_"+seed+".sh");
    elif inputMode == 'DKst0':
        os.system("echo 'bin/run -m D2PiKPi -L data/FullSel -d s21:s21r1:s24:s26 "+others+" -c DD0:DKst0 "+brstring+" -v "+variableList+" -V systematicsDir -s "+seed+"'>> BatchSubmit_"+seed+".sh");
      
    os.chdir(initDir)

def makePlotScript(sysDir, wall, brstring,others):
    initDir=os.getcwd();
    dirname=sysDir
    os.chdir(dirname)
    os.system("echo '#!/bin/bash' >> MakePlots.sh");
    os.system("echo '#PBS -l cput="+wall+":00:00' >> MakePlots.sh");
    os.system("echo '#PBS -l walltime="+wall+":00:00' >> MakePlots.sh");
    os.system("echo 'source /data/lhcb/sw/scripts/lbsetup-cvmfs.sh' >> MakePlots.sh");
    os.system("echo '. SetupProject.sh Gaudi ROOT' >> MakePlots.sh");

    os.system("echo 'cd "+dirname+"/../' >> MakePlots.sh")
    if inputMode == 'DsPhi':
        os.system("echo 'bin/run -m Ds2PhiPi:Ds2KKPi:Ds2PiPiPi:Ds2KPiPi -L data/FullSel -d s21:s21r1:s24:s26 "+others+" -c DsD0:DsPhi:DsPhiSide "+brstring+"-v  -V systematicsDir' >> MakePlots.sh");
    elif inputMode == 'DKst0':
        os.system("echo 'bin/run -m D2PiKPi -L data/FullSel -d s21:s21r1:s24:s26 "+others+" -c DD0:DKst0 "+brstring+"-v  -V systematicsDir' >> MakePlots.sh");
     
    os.chdir(initDir)

basedir = "/data/lhcb/users/hadavizadeh/B2DsPhi/Systematics"


inputMode     = str(raw_input("DsPhi or DKst0? "))
if inputMode!='DsPhi' and inputMode!= 'DKst0':
    print "Mode not recognised: " + inputMode
    print "Assuming you meant DsPhi"


desc = raw_input("Enter a description of the fits you are running (today's date will automatically be included): ")

if desc == '':
    print 'You entered a blank string so I am using "sys"'
    desc = 'sys'

sysDir = basedir + '/' + str(todaysdate) + "_" + desc + "_" + inputMode
if os.path.exists(sysDir):
    print
    print '*** WARNING: directory', sysDir, 'already exists! ***'
    decision = raw_input('*** Input Y to overwrite, anything else to abort: ').lower()
    if decision != 'y':
        sys.exit()
    else:
        shutil.rmtree(sysDir)
        print 'Directory has been removed and will be recreated'
print 'Creating', sysDir
os.mkdir(sysDir)
print


brtype = str(raw_input("Fit yields, Single branching fraction, or four branching fractions? type y, s, or f: "))
if brtype != 'y' and brtype != 'f' and brtype != 's':
    print 'Decision not recognised: '+ brtype
brstring = ''
if brtype=='s':
    brstring = ' -R '
if brtype=='f':
    brstring = ' -r '

mergedInput = str(raw_input("Merged plots? y or other: "))
if mergedInput=='y':
    brstring += ' -D '

# Get seed
seedInput = int(raw_input("Give first random seed: "))
# Get Number of batch jobs
numBatchJobs = int(raw_input("Num batch jobs: "))

print 'Possible variables:'
print '========= Signal ==========='
print '    fixedSig_n '
print '    fixedSig_alpha '
print '    fixedSig_Sigmaratio '
print '    fixedSig_Sigmafrac '
print '    fixedSig_NormSigma '
print '    fixedSig_BinRatios '
print
print '    To Double tail variations: fixedSig_double'
print 
print '    All Signal:'
print '          fixedSig_n:fixedSig_alpha:fixedSig_Sigmaratio:fixedSig_Sigmafrac:fixedSig_NormSigma:fixedSig_BinRatios' 
print
print '========= Backgrounds ==========='
print '    fixedBG_DsD0'
print '    fixedBG_DsPhi'
print '    fixedBG_Dsa1'
print '    fixedBG_hel'
print '    fixedBG_noDsstPhi'
print '    fixedBG_noDsa1'
print '    fixedBG_slope'
print '    fixedBG_Dsa1_smear'
print '    fixedBG_DsstKKst_smear'
print '    fixedBG_DD_Ratios'
print '    fixedBG_DD_Fractions'
print '    fixedBG_Comb_shape_flat'
print '    fixedBG_Comb_shape_line'
print
print '    All Background:'
print '          fixedBG_DsD0:fixedBG_DsPhi:fixedBG_Dsa1:fixedBG_hel:fixedBG_Dsa1_smear:fixedBG_DsstKKst_smear:fixedBG_slope:fixedBG_DD_Ratios:fixedBG_DD_Fractions'
print
print '========= With out some backgrounds ==========='
print '    fixedBG_noDsstPhi'
print '    fixedBG_noDsa1'
print 
print '========= Just draw plots (add to any systematic or none) ==========='
print '    draw'
print

# Get Number of Toys per job
variableList = str(raw_input("Colon seperated list of variables to vary:  "))
# Get wall time for job
numWalltime = int(raw_input("How many hours should it run for? "))

print
extraSettings = ' -S -M -Z -Q -X -H'
print "Default extra settings: " + extraSettings
changeextra = str(raw_input("Change any? type y to change, anything else to leave... "))
if changeextra.lower() == 'y':
    extraSettings = ' -M -Z -Q -X '
    changesplit = str(raw_input("Split by charge? type y to change, anything else to leave... "))
    if changesplit.lower() != 'y':
        extraSettings = extraSettings + ' -S '
    changehelsplit = str(raw_input("Split by helicity? type y to change, anything else to leave... "))
    if changehelsplit.lower() == 'y':
        extraSettings = extraSettings + ' -H '
    justDraw = str(raw_input("Just Draw pdfs? type y to change, anything else to leave... "))
    if justDraw.lower() == 'y':
        extraSettings = extraSettings + ' -0 '
    binned = str(raw_input("Binned fit? type y to change, anything else to leave... "))
    if binned.lower() == 'y':
        extraSettings = extraSettings + ' -B '

print "Using extra settings: " + extraSettings
print 
print 'Copying the code into the new toy directory.' 
# Make a copy of the code there ready to run - then don't have to wait for each toy to complete
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy')
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy/src')
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy/scripts')
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy/bin')
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy/pdfs')
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy/systematicsDir')
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy/systematicsDir/plots')
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy/batch')
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy/roodatasets')
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy/data')
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy/results')
os.system('mkdir '+sysDir+'/Bu2DsPhi_Fitter_copy/data/FullSel')
os.system("cp -r /home/hadavizadeh/Bc_Analysis/DataStripping/B_PhiD/Git_Fit/Test/DsPhiFitter/roodatasets/*.root "+sysDir+"/Bu2DsPhi_Fitter_copy/roodatasets")
os.system("cp -r /home/hadavizadeh/Bc_Analysis/DataStripping/B_PhiD/Git_Fit/Test/DsPhiFitter/src/* "+sysDir+"/Bu2DsPhi_Fitter_copy/src")
os.system("cp -r /home/hadavizadeh/Bc_Analysis/DataStripping/B_PhiD/Git_Fit/Test/DsPhiFitter/scripts/* "+sysDir+"/Bu2DsPhi_Fitter_copy/scripts")
os.system("cp -r /home/hadavizadeh/Bc_Analysis/DataStripping/B_PhiD/Git_Fit/Test/DsPhiFitter/bin/* "+sysDir+"/Bu2DsPhi_Fitter_copy/bin")
os.system("cp -r /home/hadavizadeh/Bc_Analysis/DataStripping/B_PhiD/Git_Fit/Test/DsPhiFitter/pdfs/* "+sysDir+"/Bu2DsPhi_Fitter_copy/pdfs")


jobCounter=1
sysDirName = sysDir+"/Bu2DsPhi_Fitter_copy"
os.system('cd '+sysDirName + '/batch')

for toySeed in range(seedInput,numBatchJobs+seedInput):
    makeSubmitScript(sysDirName+ '/batch', str(toySeed), variableList ,str(numWalltime), brstring, inputMode,extraSettings) 
    
    print 'Submitting job '+str(jobCounter)+" with seed "+str(toySeed)
    time.sleep(0.1)
    os.chdir(sysDirName + '/batch')
    os.system('qsub -N Systematic_'+str(toySeed)+' BatchSubmit_'+str(toySeed)+'.sh')
    print
    jobCounter+=1

makePlotScript(sysDirName+'/batch',str(numWalltime), brstring,extraSettings)
print 
print "Toys can be found in:"
os.system('pwd')
