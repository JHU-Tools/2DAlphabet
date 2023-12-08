import subprocess, json, ROOT, os, copy, time, glob
from collections import defaultdict
from contextlib import contextmanager

# Function stolen from https://stackoverflow.com/questions/9590382/forcing-python-json-module-to-work-with-ascii
def open_json(f):
    '''Open a JSON file. Specify twoDconfig to true if this is a 2DAlphabet 
    configuration file.

    Function adapted from https://stackoverflow.com/questions/9590382/forcing-python-json-module-to-work-with-ascii

    Args:
        f (str): File name and path.
        twoDconfig (bool, optional): Set to True if the JSON file is a 2DAlphabet
        configuration file. Defaults to True.

    Returns:
        dict: JSON opened as a python dictionary.
    '''
    with open(f) as fInput_config:
        input_config = json.load(fInput_config, object_hook=ascii_encode_dict)  # Converts most of the unicode to ascii

    return input_config

def ascii_encode_dict(data): 
    '''Convert a unicode encoded dictionary into ascii.

    Args:
        data (dict): Input dictionary.

    Returns:
        dict: Dict encoded with ascii instead of unicode.
    '''
    def _ascii_encode(x):
        if isinstance(x, unicode):
            return x.encode('ascii')
        elif isinstance(x, dict):
            return ascii_encode_dict(x)
        elif isinstance(x, list):
            out = []
            for s in x:
                if isinstance(s, unicode):
                    out.append(s.encode('ascii'))
                else:
                    out.append(s)
            return out
        else:
            return x

    return dict(map(_ascii_encode, pair) for pair in data.items())

def arg_dict_to_list(indict):
    '''Convert dictionary of arguments into a list of
    `<key>=<value>` strings.

    Args:
        indict (dict): Dictionary of arguments where keys are
        the argument name and values are the value for the argument to take.

    Returns:
        list(str): List of `<key>=<value>` strings
    '''
    return ['%s=%s'%(k,v) for k,v in indict.items()]

def parse_arg_dict(parser,indict):
    '''Use ArgumentParser to parse arguments in an input dictionary
    where the key is the arg name and the value is the value for the arg to take.

    Args:
        parser (ArgumentParser): Object to update with indict key-values.
        indict (dict): Dictionary holding argument names (keys) and values (values).

    Returns:
        Namespace: Object storing each argument as an attribute of the Namespace.
    '''
    new_namespace = parser.parse_args([])
    new_namespace.__dict__.update(indict)
    return new_namespace

def execute_cmd(cmd,dryrun=False,out=None): # pragma: no cover
    '''Print and execute a command-line command via subprocess.call().
    If dryrun==True, only print.

    Args:
        cmd (str): Command to execute as a subprocess.
        dryrun (bool, optional): If True, only print. Defaults to False.
    '''
    print ('Executing: '+cmd)
    if not dryrun:
        if out: cmd = cmd + ' | tee %s'%out
        subprocess.call([cmd],shell=True)

@contextmanager
def cd(newdir): # pragma: no cover
    '''Change active directory so that the local directory becomes newdir.
    This affects everything from subprocess calls to ROOT Print() and SaveAs()
    commands.

    Example:
        ::
        
            with cd('/path/to/stuff/'):
                ...
                <code acting in /path/to/stuff/ directory>
                ...

    Args:
        newdir (str): Directory to cd to within the `with` statement.
    '''
    print ('cd '+newdir)
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

def make_RDH(myTH2,RAL_vars,altname=''):
    '''Create a RooDataHist from the input TH2 and RooArgList of
    axis variables.

    Args:
        myTH2 (TH2): Histogram to turn into RooDataHist.
        RAL_vars (RooArgList): List of RooRealVars representing the axes.
        altname (str, optional): Alternate name. Defaults to '' in which case the TH2 name is used.

    Returns:
        RooDataHist
    '''
    name = myTH2.GetName() if altname == '' else altname
    thisRDH = ROOT.RooDataHist(name,name,RAL_vars,myTH2)
    return thisRDH

# def make_RHP(myRDH,RAL_vars):
#     name = myRDH.GetName()
#     thisRAS = ROOT.RooArgSet(RAL_vars)
#     thisRHP = ROOT.RooHistPdf(name,name,thisRAS,myRDH)
#     return thisRHP

def dict_copy(inDict,structureOnly=False):
    '''Recursively copy a dictionary with the option to only copy the
    dictionary structure.

    Args:
        inDict (dict): Input dictionary.
        structureOnly (bool, optional): If True, only copy the structure of the dictionary
        with stand-in value of 0. Defaults to False.

    Returns:
        dict: Copy.
    '''
    newDict = {}
    for k1,v1 in inDict.items():
        if isinstance(v1,dict):
            newDict[k1] = dict_copy(v1,structureOnly)
        else:
            if structureOnly: newDict[k1] = 0
            else: newDict[k1] = v1
    return newDict

def nested_dict(level,t):
    '''Create an empty dictionary nested to <level> levels
    and with type <t> as the default stand-in in the inner-most dictionary.

    Args:
        level (int): Number of nested levels (ie. max number of keys that can be used).
        t (type): Type to provide to defaultdict to determine the stand-in values.

    Returns:
        [type]: [description]
    '''
    if level > 1:
        out = defaultdict(lambda: nested_dict(level-1,t))
    else:
        out = defaultdict(t)
    return out

def roofit_form_to_TF1(RFVform,shift=0): # shift tells function how much to shift the indices of the coefficients by
    '''Convert a RooFit formula (using @) to a ROOT.TF1 type formula (using []).

    Args:
        RFVform (str): RooFit formula string.
        shift (int, optional): Number of indices to shift the coefficients in the case of abormal
        indexing. Defaults to 0.

    Returns:
        str: ROOT.TF1 formula string.
    '''
    TF1form = ''
    lookingForDigit = False
    for index,char in enumerate(RFVform):
        if char == '@':
            TF1form+='['
            lookingForDigit = True
        elif lookingForDigit:
            if char.isdigit():
                TF1form+=str(int(char)+shift)
                if index == len(RFVform)-1:
                    TF1form+=']'
            else:
                TF1form+=']'+char
                lookingForDigit = False
        else:
            TF1form+=char

    return TF1form

def set_hist_maximums(histList, factor=1.1):
    '''Take in a list of histograms and set the maximum of each
    to the maximum of the group multiplied by `factor`.

    Args:
        histList (list(TH1)): List of histograms to compare and set.
        factor (float, optional): Defaults to 1.1.

    Returns:
        list(TH1): Histograms with maximum set.
    '''
    # Set the max of the range so we can see all three histograms on the same plot
    out = []
    yMax = get_hist_maximum(histList)
    for h in histList:
        h.SetMaximum(yMax*factor)
        out.append(h)
    return out

def get_hist_maximum(histList):
    yMax = histList[0].GetMaximum()
    for h in range(1,len(histList)):
        if histList[h].GetMaximum() > yMax:
            yMax = histList[h].GetMaximum()
    return yMax

# def get_config_dirs(projPath):
#     '''Get the sub-directories in the project directory
#     that correspond to the configs used.

#     Args:
#         projPath (str): Path to project directory.

#     Returns:
#         list(str): List of sub-directories.
#     '''
#     return [f.path for f in os.scandir(projPath) if f.is_dir()]

def is_filled_list(d,key):
    '''Checks if the dictionary (`d`) entry at `key` is
    a non-empty list. If the key does not exist in the dictionary,
    return False. If the value in the dictionary is not a list, return False.

    Args:
        d (dict): Dictionary.
        key (non-enumerable): Key in dictionary.

    Returns:
        bool: True if `d[key]` is a non-empty list. Otherwise, False.
    '''
    if not isinstance(d, dict):
        raise TypeError('Arg d is not a dictionary.')
    return (key in d and isinstance(d[key],list) and len(d[key]) > 0)

def copy_update_dict(d1,d2):
    out = copy.deepcopy(d1)
    out.update(d2)
    return out

def replace_multi(s,findreplace):
    for f,r in findreplace.items():
        if f in s:
            s = s.replace(f,r)
    return s

def unpack_to_line(toUnpack):
    return ' '.join(['{%s:20}'%i for i in range(len(toUnpack))]).format(*toUnpack)

def _combineTool_impacts_fix(fileNameExpected):
    # Ex. higgsCombine_initialFit_Test.MultiDimFit.mH0.root is needed but higgsCombine_initialFit_Test.MultiDimFit.mH0.123456.root is created if a toy is used.
    seed_version = '%s.*.root'%('.'.join(fileNameExpected.split('.')[:-1]))
    potential_files_to_rename = glob.glob(seed_version)

    # Only run if there are seeded files to rename
    if potential_files_to_rename > 0:
        all_seeds = list(set([f.split('.')[-2] for f in potential_files_to_rename]))
        if len(all_seeds) > 1:
            raise RuntimeError('More than one seed found when trying to move files for combineTool (%s). Clean up the area and try again.'%all_seeds)
        
        for f in potential_files_to_rename:
            new_name = '%s.root'%('.'.join(f.split('.')[:-2]))
            execute_cmd('mv %s %s'%(f,new_name))

# ----------------- Inline condor submission --------------------
class CondorRunner():
    def __init__(self, name, primaryCmds, toPkg, runIn, toGrab, remakeEnv=False, eosRootfileTarball=None):
        '''Should be run in a CMSSW/src directory and not in a nested folder. All paths
        (besides the EOS rootfile tarball path) are treated relative to this folder.

        Args:
            name ([type]): [description]
            primaryCmds ([type]): [description]
            toPkg ([type]): [description]
            runIn ([type]): [description]
            toGrab ([type]): [description]
            remakeEnv (bool, optional): [description]. Defaults to False.
            eosRootfileTarball ([type], optional): [description]. Defaults to None.
        '''
        self.name = name
        self.run_in = runIn
        self.to_grab = toGrab
        self.primary_cmds = primaryCmds
        self.rootfile_tarball_path = eosRootfileTarball
        self.cmssw = os.environ['CMSSW_BASE'].split('/')[-1]
        if not os.path.exists('notneeded/'): execute_cmd('mkdir notneeded')
        
        self.env_tarball_path = make_env_tarball(remakeEnv)
        self.pkg_tarball_path = self._make_pkg_tarball(toPkg) if toPkg != None else ''
        self.run_script_path = self._make_run_script()
        self.run_args_path = self.run_script_path.replace('.sh','_args.txt')

    def submit(self):
        abs_path = os.getcwd()
        abs_twoD_dir_base = abs_path[:abs_path.find(self.cmssw)]+self.cmssw+'/src/2DAlphabet/TwoDAlphabet/'
        timestr = time.strftime("%Y%m%d-%H%M%S")
        out_jdl = 'temp_'+timestr+'_jdl'

        execute_cmd("sed 's$TEMPSCRIPT${0}$g' {1}/condor/jdl_template > {2}".format(self.run_script_path, abs_twoD_dir_base, out_jdl))
        execute_cmd("sed -i 's$TEMPTAR${0}$g' {1}".format(self.pkg_tarball_path, out_jdl))
        execute_cmd("sed -i 's$TEMPARGS${0}$g' {1}".format(self.run_args_path, out_jdl))
        execute_cmd("condor_submit "+out_jdl)
        execute_cmd("mv {0} notneeded/".format(out_jdl))

    def _make_pkg_tarball(self,to_pkg):
        start_dir = os.getcwd()
        out_dir = start_dir+'/notneeded'
        out_path = '%s/%s_input.tgz'%(out_dir,self.name)
        with cd(os.environ['CMSSW_BASE']+'/src'):
            if os.path.exists(out_path):
                execute_cmd('rm '+out_path)
            print ('Making package tarball %s.tgz'%self.name)
            execute_cmd('tar --exclude=*.tgz -czf {0}.tgz {1}'.format(self.name, to_pkg))
            print ('Done')
            execute_cmd('mv %s.tgz %s'%(self.name,out_path))

        return out_path

    def _make_run_script(self):
        blocks = []
        blocks.append(
            _setup_env.format(
                env_tarball_path=self.env_tarball_path,
                cmssw=self.cmssw,
                scram_arch=os.environ['SCRAM_ARCH']
            )
        )

        if self.pkg_tarball_path != None:
            blocks.append(_setup_tar_pkg.format(pkg_tarball=self.pkg_tarball_path.split('/')[-1], cmssw=self.cmssw))
        
        blocks.append('cd %s/src/; eval `scramv1 runtime -sh`'%self.cmssw)

        if self.rootfile_tarball_path != None:
            blocks.append(_unpack_rootfiles.format(rootfile_tarball_path=self.rootfile_tarball_path))

        blocks.append('cd %s'%self.run_in)
        blocks.append('echo $*')
        blocks.append('$*')
        blocks.append('cd $CMSSW_BASE/src/')
        blocks.append(_grab_output.format(out_id='%s_output_${CONDOR_ID}'%(self.name), to_grab=self.to_grab))

        shell_name = 'run_'+self.name+'.sh'
        with open(shell_name,'w') as run_script:
            for block in blocks:
                run_script.write(block+'\n')

        with open(shell_name.replace('.sh','_args.txt'),'w') as run_args:
            for c in self.primary_cmds:
                run_args.write(c+'\n')

        return os.path.abspath(shell_name)

def make_env_tarball(makeEnv=True):
    dir_base = os.environ['CMSSW_BASE']
    cmssw = dir_base.split('/')[-1]
    user = os.environ['USER']
    out_eos_path = 'root://cmseos.fnal.gov//store/user/{user}/{cmssw}_env.tgz'.format(user=user,cmssw=cmssw)
    if makeEnv:
        with cd(dir_base+'/../'):
            if os.path.exists('%s_env.tgz'%cmssw):
                execute_cmd('rm %s_env.tgz'%cmssw)
            print ('Making env tarball %s_env.tgz...'%cmssw)
            execute_cmd('tar --exclude-caches-all --exclude-vcs --exclude-caches-all --exclude-vcs -czf {cmssw}_env.tgz {cmssw} --exclude=tmp --exclude=".scram" --exclude=".SCRAM"'.format(cmssw=cmssw))
            print ('Done')
            execute_cmd('xrdcp {cmssw}_env.tgz {out}'.format(cmssw=cmssw,out=out_eos_path))
    
    return out_eos_path

# Done in base dir
_setup_env = '''#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
xrdcp {env_tarball_path} env_tarball.tgz
export SCRAM_ARCH={scram_arch}
scramv1 project CMSSW {cmssw}
tar -xzf env_tarball.tgz
rm env_tarball.tgz'''

# Done in base dir
_setup_tar_pkg = '''
mkdir tardir; cp {pkg_tarball} tardir/; cd tardir
tar -xzvf {pkg_tarball}
rm {pkg_tarball}
cp -r * ../{cmssw}/src/
cd ../'''

# Done in CMSSW/src
_unpack_rootfiles = '''
mkdir EOS_ROOTFILES
cd EOS_ROOTFILES
xrdcp {rootfile_tarball_path} ./
tar -xzf {rootfile_tarball}
rm {rootfile_tarball}
cd ../'''

# Done in CMSSW/src
_grab_output = '''
tar -czvf {out_id}.tgz {to_grab}
cp {out_id}.tgz $CMSSW_BASE/../'''
