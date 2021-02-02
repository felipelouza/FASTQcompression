#!/usr/bin/env python3

import sys, time, argparse, subprocess, os.path, shutil, struct

Version = "v0.1"


Description = """Tool to compress FASTQ files

For example
  {exe} INPUT.fastq -o OUTPUT -1
will produce the files OUTPUT.fq.7z
 
--------------------------
Command line options:
--------------------------
""".format(exe=sys.argv[0])

gsufsort_exe = "external/gsufsort/gsufsort"
header_split = "sed -n 1~4p"
qs_split = "sed -n 4~4p"
dna_split = "sed -n 2~4p"
gzip_exe = "gzip -9 -k -f"
#zip7_exe = "7z a -mx9 -mmt12"
zip7_exe = "7z a -mm=PPMd"

smooth_exe = "src/fq_compression"

def main():
    parser = argparse.ArgumentParser(description=Description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('input',      help='input file name', type=str, nargs='+')
    parser.add_argument('-o','--out', help='output base name (def. input base name)', default="", type=str)  
    #parser.add_argument('--delete', help='delete output files',action='store_true')
    parser.add_argument('--step1',    help='stop after step 1 (debug only)',action='store_true')
    parser.add_argument('--step2',    help='stop after step 2 (debug only)',action='store_true')
    parser.add_argument('--step3',    help='stop after step 3 (debug only)',action='store_true')
    parser.add_argument('--step4',    help='stop after step 4 (debug only)',action='store_true')
    parser.add_argument('--original', help='do not call step 2',action='store_true')
    parser.add_argument('-1', '--m1', help='mode 1: FASTQ', action='store_true',default=True)
    parser.add_argument('-2', '--m2', help='mode 2: BWT+QS', action='store_true')
    parser.add_argument('-3', '--m3', help='mode 3: BWT+QS+H', action='store_true')
    parser.add_argument('--no_headers',  help='ignores the headers', action='store_true', default=False)
    parser.add_argument('--others',   help='run all competitors',action='store_true')
    parser.add_argument('-v',         help='verbose: extra info in the log file',action='store_true')
    args = parser.parse_args()
    # ---- check number of input files and define basename
    check_input(args)
    # ---- create and open log file
    logfile_name = args.basename + ".log"
    # get main directory
    args.dir = os.path.split(sys.argv[0])[0]
    print("Sending logging messages to file:", logfile_name)

    with open(logfile_name,"w") as logfile:
        
        ##
        if(args.m2): args.m1 = args.m3 = False
        if(args.m3): args.m1 = args.m2 = False
        ##
        if(not args.m1): args.no_headers = True
        ##

        print(">>> fastq-bwt version " + Version,file=logfile)
        print(">>> fastq-bwt version " + Version)
        if(args.m1):
            print(">>> mode 1: FASTQ",file=logfile) 
            print(">>> mode 1: FASTQ") 
        if(args.m2):
            print(">>> mode 2: BWT+QS",file=logfile) 
            print(">>> mode 2: BWT+QS")
        if(args.m3):
            print(">>> mode 3: BWT+QS+H",file=logfile) 
            print(">>> mode 3: BWT+QS+H")

        show_command_line(logfile)
        logfile.flush()

        if len(args.out)==0 : args.out=args.input[0]

        # temporary files
        args.tmp = []
        args.stream = []

        #--- step1: compute BWT+QS
        start = time.time()
        if(step1(args, logfile, logfile_name)!=True):
            sys.exit(1)
        print("Elapsed time: {0:.4f}".format(time.time()-start))
        if args.step1:
            print("Exiting after step 1 as requested")
            return
    
        # --- step2: smooth BWT and QS sequences 
        start = time.time()
        if(step2(args, logfile, logfile_name)!=True):
            sys.exit(1)
        print("Elapsed time: {0:.4f}".format(time.time()-start))
        if args.step2:
            print("Exiting after step 2 as requested")
            return

        if(args.m3):
            # --- step3: extract headers
            start = time.time()
            if(step3(args, logfile, logfile_name)!=True):
                sys.exit(1)
            print("Elapsed time: {0:.4f}".format(time.time()-start))
            if args.step3:
                print("Exiting after step 3 as requested")
                return

        # --- step4: compute BWT+QS
        if(args.m2 or args.m3):
            start = time.time()
            if(step4(args, logfile, logfile_name)!=True):
                sys.exit(1)
            print("Elapsed time: {0:.4f}".format(time.time()-start))
            if args.step1:
                print("Exiting after step 4 as requested")
                return
        
        # compressed files
        args.output = []

        # --- step5: compress new FASTQ, new ebwt string, new qs string and headers separatedly 
        start = time.time()
        if(step5(args, logfile, logfile_name)!=True):
            sys.exit(1)
        print("Elapsed time: {0:.4f}".format(time.time()-start))

        # ---- final report
        insize = os.path.getsize(args.input[0])

        print("=== results ==="); 
        print("Original:\t{0:.2f} MB".format(insize/(1024*1024)))
        if(args.v): print(args.input[0])
        print("==")
        outsize = 0
        for f in args.output:
            outsize += os.path.getsize(f)
        print("Compressed:\t{0:.2f} MB".format(outsize/(1024*1024)))
        print("Ratio = {0:.2f}".format(outsize/insize))
        if(args.v):
            for f in args.output:
               print(f) 

        # --- extra: compress INPUT.fastq with default method 
        if not args.others:
            return

        #gzip
        start = time.time()
        if(gzip(args, logfile, logfile_name)!=True):
            sys.exit(1)
        print("Elapsed time: {0:.4f}".format(time.time()-start))
        outsize = os.path.getsize(args.input[0]+".gz")
        print("Compressed:\t{0:.2f} MB".format(outsize/(1024*1024)))
        print("Ratio = {0:.2f}".format(outsize/insize))


        #7zip
        start = time.time()
        if(zip7(args, logfile, logfile_name)!=True):
            sys.exit(1)
        print("Elapsed time: {0:.4f}".format(time.time()-start))
        outsize = os.path.getsize(args.input[0]+".7z")
        print("Compressed:\t{0:.2f} MB".format(outsize/(1024*1024)))
        print("Ratio = {0:.2f}".format(outsize/insize))



    return True

##

def step1(args, logfile, logfile_name):
    print("--- Step 1 ---", file=logfile); logfile.flush()
    exe = os.path.join(args.dir, gsufsort_exe)
    options = ""
    if len(args.out)>0 : options+="-o "+args.out
    else : options+="-o "+args.input[0]
    command = "{exe} {ifile} --bwt --qs {opt}".format(exe=exe, ifile=args.input[0], opt=options)
    print("=== gsufsort ==="); print(command)
    # tmp files
    args.tmp.append(args.basename+".bwt")
    args.tmp.append(args.basename+".bwt.qs")
    return execute_command(command, logfile, logfile_name)

def step2(args, logfile, logfile_name):
    #TODO: check
    if args.original:
        print("--- Step 2 ---", file=logfile); logfile.flush()
        command = "cp "+ args.input[0] +" "+args.out+".fq" 
        print(command)
        os.system(command)
    else:
        print("--- Step 2 ---", file=logfile); logfile.flush()
        exe = os.path.join(args.dir, smooth_exe)
        options = "-e " + args.tmp[0] + " -q " + args.tmp[1] + " -f " + args.input[0]+" -o "+args.out+".fq"
        if(args.no_headers): #ignore headers
            options+=" -H"
        command = "{exe} {opt}".format(exe=exe, opt=options)
        print("=== smooth-qs ===")
        print(command)
        if(args.m1): args.stream.append(args.out+".fq")
        return execute_command(command, logfile, logfile_name)
    return True

##
def step3(args, logfile, logfile_name):
    print("--- Step 3 ---", file=logfile); logfile.flush()
    ##
    exe = header_split
    ifile = args.input[0]
    ofile = args.out+".h"
    command = "{exe} {ifile} > {ofile}".format(exe=exe, ifile=ifile, ofile=ofile)
    print("=== header ===")
    print(command)
    os.system(command)
    args.stream.append(args.basename+".h")
    return True
##

def step4(args, logfile, logfile_name):
    print("--- Step 4 ---", file=logfile); logfile.flush()
    exe = os.path.join(args.dir, gsufsort_exe)
    options = ""
    if len(args.out)>0 : options+="-o "+args.out+".fq"
    else : options+="-o "+args.out+".fq"
    command = "{exe} {ifile} --bwt --qs {opt}".format(exe=exe, ifile=args.out+".fq", opt=options)
    print("=== gsufsort ==="); print(command)
    # stream files
    args.stream.append(args.out+".fq.bwt")
    args.stream.append(args.out+".fq.bwt.qs")
    return execute_command(command, logfile, logfile_name)

def step5(args, logfile, logfile_name):
    print("--- Step 5 ---", file=logfile); logfile.flush()
    #exe = gzip_exe
    exe = zip7_exe
    print("=== compression ===")
    for f in args.stream:
        #ofile = f+".gz"
        #command = "{exe} {ifile}".format(exe=exe, ifile=f)
        ofile = f+".7z"
        command = "{exe} {ofile} {ifile}".format(exe=exe, ifile=f, ofile=ofile)
        print(command)
        execute_command(command, logfile, logfile_name)
        args.output.append(ofile)
    return True

def gzip(args, logfile, logfile_name):
    print("--- gzip ---", file=logfile); logfile.flush()
    exe = gzip_exe
    print("=== gzip ===")
    command = "{exe} {ifile}".format(exe=exe, ifile=args.input[0])
    print(command)
    return execute_command(command, logfile, logfile_name)

def zip7(args, logfile, logfile_name):
    print("--- 7z ---", file=logfile); logfile.flush()
    exe = zip7_exe
    ofile = args.input[0]+".7z"
    print("=== 7z ===")
    command = "{exe} {ofile} {ifile}".format(exe=exe, ifile=args.input[0], ofile=ofile)
    print(command)
    return execute_command(command, logfile, logfile_name)

########

# check correctness of number of input file and define basename for output
def check_input(args):
    if len(args.out)==0:       # specify basename for input files gap+merge
        args.basename = args.input[0]
    else:
        args.basename = args.out
    return True

# compute hash digest for a file 
def file_digest(name,logfile):
    try:
        hash_command = "{exe} {infile}".format(exe=shasum_exe, infile=name)
        hashsum = subprocess.check_output(hash_command.split(),stderr=logfile)
        hashsum = hashsum.decode("utf-8").split()[0]
    except:
        hashsum = "Error!" 
    return hashsum  

# execute command: return True is everything OK, False otherwise
def execute_command(command,logfile,logfile_name):
    try:
        subprocess.check_call(command.split(),stdout=logfile,stderr=logfile)
    except subprocess.CalledProcessError:
        print("Error executing command line:")
        print("\t"+ command)
        print("Check log file: " + logfile_name)
        return False
    return True

def show_command_line(f):
    f.write("Python command line: ") 
    for x in sys.argv:
        f.write(x+" ")
    f.write("\n")   

if __name__ == '__main__':
    main()
