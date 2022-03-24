__author__ = 'djk'

from tkinter import *
import sys
import os
import subprocess
from multiprocessing.pool import ThreadPool as Pool
import shutil


# TSC1 9:135766735-135820020
# TSC2 16:2097990-2138713
# MTOR 1:11166588-11322608
# PTEN: 10:89623195-89728532
# PIK3CA: 3:178866311-178952497
# VHL: 3:10183319-10195354
# FLCN 17:17113000-17143000
# TXNIP: "1:145438462-145442628"
# KDM6A: "X:44732423-44971845"
# MSR1: "8:15965387-16050300"
# NF2: "22:29999545-30094589"
# VHL: "3:10183319-10195354"
# NPNT: "4:106816597-106892828"
# MET: "7:116312459-116438440"
# TP53: "17:7571720-7590868"
# BAP1: "3:52435020-52444121"
# PTEN: "10:89623195-89728532"
# ATM: "11:108093559-108239826"
# MICALCL: "11:12308447-12380691"
# CDKN1A: "6:36644237-36655116"
# HIF1A: "14:62162119-62214977"
# TSC1: "9:135766735-135820020"
# SMARCB1: "22:24129150-24176705"
# TSC2: "16:2097990-2138713"
# RHEB: "7:151163098-151217010"
# PIK3CA: "3:178866311-178952497"
# PBRM: "3:52579368-52713739"
# TCEB1: "8:74857373-74884522"
# MTOR: "1:11166588-11322608"
# NFE2L2: "2:178095031-178129859"
# SLITRK6: "13:86366922-86373483"
# SETD2: "3:47057898-47205467"
# STAG2: "X:123094475-123236505"
# KDM5C: "X:53220503-53254604"

class Application(Frame):
 
    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.pack()
        self.createWidgets()

    def createWidgets(self):
        '''
        Creates the buttons and entry boxes for the user input for the analysis
        '''
        self.genbai = Button(self, text="Generate .bai", command=self.baigenerator)
        self.genbai.pack(side=TOP)

        L1 = Label(self, text="Location (chr:nt-nt):")
        L1.pack()
        global E1
        # E1 = Entry(self, bd =5)
        # E1.pack(side=LEFT)

        # Change the entry to a drop-down widget from which we can select a location in a specified chromosome
        E1 = StringVar(self)

        E1.set("9:135766735-135820020")
        option = OptionMenu(self, E1, "1:145438462-145442628", \
                                        "X:44732423-44971845", \
                                        "8:15965387-16050300", \
                                        "22:29999545-30094589", \
                                        "3:10183319-10195354", \
                                        "4:106816597-106892828", \
                                        "7:116312459-116438440", \
                                        "chr17:7571720-7590868", \
                                        "3:52435020-52444121", \
                                        "10:89623195-89728532", \
                                        "11:108093559-108239826", \
                                        "11:12308447-12380691", \
                                        "6:36644237-36655116", \
                                        "14:62162119-62214977", \
                                        "9:135766735-135820020", \
                                        "22:24129150-24176705", \
                                        "chr16:2097990-2138713", \
                                        "7:151163098-151217010", \
                                        "3:178866311-178952497", \
                                        "3:52579368-52713739", \
                                        "8:74857373-74884522", \
                                        "1:11166588-11322608", \
                                        "2:178095031-178129859", \
                                        "13:86366922-86373483", \
                                        "3:47057898-47205467", \
                                        "X:123094475-123236505", \
                                        "X:53220503-53254604")

        option.pack()

        L2 = Label(self, text="Minimum Variant Allele Frequency:")
        L2.pack()
        global E2
        E2 = Entry(self, bd=5)
        E2.pack()

        L3 = Label(self, text="Minimum read count (total):")
        L3.pack()
        global E3
        E3 = Entry(self, bd=5)
        E3.pack()

        L4 = Label(self, text="Minimum indel frequency:")
        L4.pack()
        global E4
        E4 = Entry(self, bd=5)
        E4.pack()

        self.start = Button(self, text="Begin Analysis", command=self.combinefunctions)
        self.start.pack()

        self.matlab = Button(self, text="Run matlab script", command=self.runmatlab)
        self.matlab.pack()

        # self.spreadsheet = Button(self, text="Excel Spreadsheet", command=self.createspreadsheet)

    def saveminfreq(self):
        '''
        Saves the input value for the minimum variant allele frequency into a text file, which will be imported as
        a variable in matlab.
        '''

        minfreq = E2.get()
        minfreqfile = open("minfreq.txt", "w")
        minfreqfile.write("%0.5f" % float(minfreq))
        minfreqfile.close()

    def saveminreadcount(self):
        '''
        Saves the input value for the minimum read count into a text file, which will be imported as a variable
        in matlab.
        '''

        minreadcount = E3.get()
        minreadcountfile = open("minreadcount.txt", "w")
        minreadcountfile.write("%d" % int(minreadcount))
        minreadcountfile.close()

    def saveindelfreq(self):
        '''
        Saves the input value for the minimum indel frequency into a text file, which will be imported as a variable
        in matlab.
        '''
        minIndelFreq = E4.get()
        minIndelFreqFile = open("minIndelFreq.txt", "w")
        minIndelFreqFile.write("%0.5f" % float(minIndelFreq))
        minIndelFreqFile.close()

    def runmatlab(self):
        '''
        Runs Matlab_commands_3_16_GenomADupdate.m script
        '''
        subprocess.call(['matlab -nojvm -r Matlab_commands_3_16_GenomADupdate'], shell=True)
        print('Matlab script completed. Open the file aafilterdata.csv\n')
        exit()

    def combinefunctions(self):
        '''
        Combines the functions that generate and modify the files in to a
        single function which can be invoked when the "Begin Analysis" button
        is pressed
        '''

        self.saveminfreq()
        self.saveminreadcount()
        self.saveindelfreq()

        self.namereader()

        if (not self.isvalidregion(geneloc)):
            print('Error, region is not valid\n')
            exit()
        else:
            print('Valid region.\n')

        self.pupgenerator()

        print("Pileup files generated.\n")

        self.outgenerator()

        print("Conversion output files generated.\n")

        self.agenerator()

        print("bam.a files generated.\n")

        self.zeroscreator()

        # Import data from GnomAD Browser
        source = open('genomeADalts.txt', 'r')
        target = open('reformatgenomeAD.txt', 'w')

        for line in source.readlines():
            if (line == 'A\n'):
                target.write('1\n')
            elif (line == 'G\n'):
                target.write('2\n')
            elif (line == 'C\n'):
                target.write('3\n')
            elif (line == 'T\n'):
                target.write('4\n')
            else:
                for char in line:
                    if (char == 'A'):
                        target.write('1')
                    elif (char == 'G'):
                        target.write('2')
                    elif (char == 'C'):
                        target.write('3')
                    elif (char == 'T'):
                        target.write('4')
                target.write('\n')

        source.close()
        target.close()

        print("Python analysis complete. mergez.txt file and namelist.txt file generated.\n")

    def isvalidregion(self, geneloc):
        '''
        Compares the user input region values with the indices of the human genome reference values and makes sure
        that the start and end index are plus or minus 10000 of the reference
        '''
        i = geneloc.find(':') + 1
        j = geneloc.find('-')

        inputLen = len(geneloc)

        inputStart = int(geneloc[i: j])

        inputEnd = int(geneloc[j + 1: inputLen])

        refSeqFile = open("genomewithinst.txt", "r")

        fileData = refSeqFile.readline()

        indexbegin = fileData.find(':') + 1
        indexmid = fileData.find('-')
        indexend = fileData.find('\'') - 2

        refStart = int(fileData[indexbegin: indexmid])
        refEnd = int(fileData[indexmid + 1: indexend])

        text_file = open("initindex.txt", "w")
        text_file.write("%d" % refStart)
        text_file.close()

        text_file2 = open("endindex.txt", "w")
        text_file2.write("%d" % refEnd)
        text_file2.close()

        target_file = open("genome.txt", "w")
        shutil.copyfileobj(refSeqFile, target_file)

        refSeqFile.close()
        target_file.close()

        openFile = open("genome.txt", "r")
        targetFile = open("genomenums.txt","w")

        for line in openFile.readlines():
            for char in line:
                if(char == "A"):
                    targetFile.write("1\n")
                elif(char == "G"):
                    targetFile.write("2\n")
                elif(char == "C"):
                    targetFile.write("3\n")
                elif(char == "T"):
                    targetFile.write("4\n")

        openFile.close()
        targetFile.close()

        return ((refStart + 10000 == inputStart) and (refEnd - 10000 == inputEnd))

    def namereader(self):
        '''
        Reads through a list of files to get the bam files and stores the
        names of these files in namelist.txt
        '''
        global geneloc
        geneloc = E1.get()

        numSamples = 0

        oldstdout = sys.stdout
        f = open('namelist.txt', 'w')
        sys.stdout = f

        # path = "/Users/guest/Desktop/Karthik"

        path = os.path.realpath(__file__)

        path = path.replace('/Analyzer.py', '')

        dirs = os.listdir(path)

        for file in dirs:
            if file.endswith(".bam"):
                # file=file.replace(".pup","")
                print(file)
                numSamples += 1

        sys.stdout = oldstdout

        numSamplesFile = open("numSamples.txt", "w")
        numSamplesFile.write("%d" % numSamples)
        numSamplesFile.close()

    def baigenerator(self):
        '''
        Generates .bai files by calling a subprocess which calls the samtools
        index command on the UNIX command line
        '''
        for line in open('namelist.txt', 'r'):
            line2 = line.replace("\n", "")
            line3 = line2.replace(" ", line2)
            subprocess.call(['samtools index ' + line3], shell=True)

    def pupgenerator(self):
        '''
        Generates .pup files by calling a subprocess which calls the samtools
        mpileup command on the UNIX command line
        '''
        for line in open('namelist.txt', 'r'):
            line2 = line.replace("\n", "")
            line3 = line2.replace(" ", line2)
            subprocess.call(['samtools mpileup -l MHPA_TSC2.bed -d0 -x -A -B ' + line3 + ' > ' + line3 + '.pup'], shell=True)

    def outgenerator(self):
        '''
        Generates .out files by calling a subprocess on the UNIX command line
        that uses the Python script v12-q50.py
        '''

        lineList = open('namelist.txt', 'r')

        numSamplesFile = open("numSamples.txt", "r")
        numSamples = int(numSamplesFile.read())

        # Create the pool for the threads
        pool = Pool(numSamples)

        numSamplesFile.close()

        def out_loop_operation(line):
            '''
            Function that represents the operation done in the loop in the
            outgenerator function
            '''
            line2 = line.replace("\n", "")
            line3 = line2.replace(" ", line2)
            # print("operation before call: %s\n" % line3)
            subprocess.call(['python v12-q50.py ' + line3 + '.pup ' + line3 + '.out'], shell=True)
            # print("operation after call: %s\n" % line3)

        # Make the threads run in parallel
        for line in lineList:
            # print("thread being run is %s\n" % line)
            pool.apply_async(out_loop_operation, (line,))

        pool.close()
        pool.join()

    def agenerator(self):
        '''
        Generates .a files by calling a subprocess which uses the cut command
        on the UNIX command line
        '''

        for line in open('namelist.txt', 'r'):
            line2 = line.replace("\n", "")
            line3 = line2.replace(" ", line2)
            subprocess.call(['cut -f 2,3,6-9,11-14,16-18 ' + line3 + '.out > ' + line3 + '.a'], shell=True)

    def zeroscreator(self):
        '''
        Populates the mergez.txt file initially with zeroes
        '''
        # rowcount = str(subprocess.check_output(['wc -l < '+largestname], shell=True))

        subprocess.call(['python zeroscreator.py zeros.txt 70000 13'], shell=True)

        for line in open('namelist.txt', 'r'):
            line2 = line.replace("\n", "")
            line3 = line2.replace(" ", line2)
            subprocess.call(['cat ' + line3 + '.a zeros.txt > ' + line3 + '.b'], shell=True)

        filenames = ''
        for line in open('namelist.txt', 'r'):
            line2 = line.replace("\n", "")
            line3 = line2.replace(" ", line2)
            filenames = filenames + line3 + '.b '
        subprocess.call(['python merger.py ' + filenames + '>>mergez.txt'], shell=True)


root = Tk()
app = Application(master=root)

app.master.title("Analyzer")

app.mainloop()
