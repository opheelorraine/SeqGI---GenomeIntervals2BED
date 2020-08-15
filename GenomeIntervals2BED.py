import sys, os, datetime
from time import strftime
from optparse import OptionParser

DefaultCOLMAP = { 		
    	'chromosome' : 0,
    	'start'      : 1,
    	'end'        : 2,
        'strand'     : None,
        'ID'         : None
	}

class GenomicInterval:

    def __init__(self,  **kwargs):
        self.fname = kwargs['fname']
        self.oname = kwargs['oname'] #os.path.splitext(self.filename)[0] + '.bed'
        
        self.SEPFILE = kwargs['sep']
        self.wstrand = kwargs['wstrand'] #Define window on strand. (Same as BEDtools "-sw")
        self.window = kwargs['window']
        self.columns = kwargs['columns']
        self.cID = kwargs['cID']
        self.onebased = kwargs['onebased']
        
        if self.columns == 'None': self.columns = None
        if self.cID == 'None': self.cID = None
        if self.window == 'None': self.window = None

        self.warnings = []
        
        self.checkargs()

    #-------------------------------------
    def checkargs(self):
        if self.cID != None:
            try: self.cID = int(self.cID)
            except ValueError:
                print("ERROR: column position (ID) has to be an integer")
                sys.exit()
                
        sepdic = {"Space": ' ', "Tab": "\t", "Comma": ",", "Semicomma": ";"}
        try: self.SEPFILE = sepdic[self.SEPFILE]
        except: pass

        self.columnMap = self.mkColumnMap(self.columns, self.cID)

        self.window  = self.readstringwindows(self.window)

        if self.window != None:
            if self.window[0] in ('TSS','TES','GeneBody') and not self.wstrand:
                warning = 'WARNING: Defining window %s relative to TSS, TES or GeneBody has to be done based on strands. "wstrand=True" was used' %self.window
                self.warnings.append(warning)
                self.wstrand = True
            
            if self.columnMap['strand'] == None and self.wstrand:
                raise NameError('ERROR: --columns=%s. Strand is needed to define windows based on strand' %self.columns)
                sys.exit()

        if self.window == None:
            warning = 'WARNING: unused argument: "wstrand=True" for no specified window.'
            self.wstrand = False
        
    #-------------------------------------
    def mkColumnMap(self, colsArg, IDcol):
        columnMap = {}
        if colsArg == None: return DefaultCOLMAP

        try:
            cols = [int(x) for x in colsArg.split(',')]

        except:
            raise NameError('ERROR: column arguments %s have to be integers' %colsArg)
            sys.exit()
            
        if len(cols) < 3:
            raise NameError('ERROR: column arguments %s have to be at least 3: chromosome,start,end' %colsArg)
            sys.exit()
            
        columnMap['chromosome'] = cols[0] - 1
        columnMap['start']      = cols[1] - 1
        columnMap['end']        = cols[2] - 1                
        if len(cols) == 4: columnMap['strand'] = cols[3] - 1
        else: columnMap['strand'] = None 
        if IDcol != None: columnMap['ID'] = IDcol - 1
        else: columnMap['ID'] = None
        
        return columnMap

    #-------------------------------------
    def readstringwindows(self,window):
        
        if window != None:
            w = window.split(':')
             
            if w[0] == "GeneBody" and len(w) < 3:
                warning = 'WARNING: "GeneBody" is normally used to define a window upstream and downstream a gene, e.g. GeneBody:-1000:1000. ' + \
                          'So, In this case the window chosen "GeneBody" was set to "None"'
                self.warnings.append(warning)
                return None #window = None
            
            try:
                w[1] = int(w[1])
                w[2] = int(w[2])
            except:
                raise NameError('ERROR: window %s not valid: window coordinates have to be integers' %w)
                sys.exit()
                
            if w[0] not in ('TSS','TES', 'GeneBody'):
                try: w[0] = int(w[0]) - 1
                except:
                    raise NameError('ERROR: window not valid: window has to be defined relative to TSS,TES,GeneBody or an integer. The integer indicates the column position where the coordinate of interest can be found in the input file')
                    sys.exit() 

            if len(w) < 3:
                raise NameError('ERROR: window %s is not valid: have to include upstream and downstream of window-relative position' %w)
                sys.exit()
                                    
            return w
        
        else: return None
        
    #-------------------------------------
    def go(self):
        '''Converts genomic coordinates to the ones specified in window'''

        self.starttime = datetime.datetime.now()
        datenow = strftime("%Y-%m-%d %H:%M:%S")
        yield 'GenomeIntervals2BED6.py'
        yield 'Message: Starting SeqGI GenomeIntervals2BED6 module'
        sys.stdout.flush()
        
        inputf = open(self.fname, 'r')
        outputf = open(self.oname, 'w')
        
        for line in inputf:
            if not (line=='\n' or line == '\r\n' or line == '\r' or line[0] == '#' or line.startswith('track') or line.startswith('browser')):
                line = line.strip().split(self.SEPFILE)
                chromosome = line[self.columnMap['chromosome']]
                try: strand = line[self.columnMap['strand']]
                except: strand = '*'
                try: ID = line[self.columnMap['ID']]
                except: ID = 'None'
                passit = False

                left = int(line[self.columnMap['start']])
                right   = int(line[self.columnMap['end']])
                if self.onebased: left = left-1

                if self.window == None:
                    left_coord = left
                    right_coord = right

                elif self.window[0] not in ('TSS', 'TES', 'GeneBody'):
                    position = int(line[self.window[0]])
                    if self.onebased: position = position-1
                    if self.wstrand:
                        if strand == '-':
                            left_coord = position + (int(self.window[2])*-1)
                            right_coord = position + (int(self.window[1])*-1)
                                
                        elif strand == '+': 
                            left_coord = position + (int(self.window[1]))
                            right_coord = position + (int(self.window[2]))
                            
                        else: passit = True

                    else:
                        left_coord = position + int(self.window[1])
                        right_coord = position + int(self.window[2])
                        
                elif self.window[0] in ('TSS', 'TES', 'GeneBody'):
                    
                    if strand == '+': (tss,tes) = (left,right)
                    elif strand == '-': (tss,tes) = (right,left)
                    else: passit = True #if there is no strand then it cannot compute interval...

                    if strand == '-':
                        if self.window[0] == 'GeneBody':
                            right_coord = tss + (int(self.window[1])*-1)
                            left_coord = tes + (int(self.window[2])*-1)
                        elif self.window[0] == 'TSS':
                            left_coord = tss + (int(self.window[2])*-1)
                            right_coord = tss + (int(self.window[1])*-1)
                        elif self.window[0] == 'TES':
                            left_coord = tes + (int(self.window[2])*-1)
                            right_coord = tes + (int(self.window[1])*-1)

                    elif strand == '+':
                        if self.window[0] == 'GeneBody':
                            left_coord = tss + (int(self.window[1]))
                            right_coord = tes + (int(self.window[2]))
                        elif self.window[0] == 'TSS':
                            left_coord = tss + int(self.window[1])
                            right_coord = tss + int(self.window[2])
                        elif self.window[0] == 'TES':
                            left_coord = tes + int(self.window[1])
                            right_coord = tes + int(self.window[2])
                    else: passit = True
                    
                        
                if not passit:            
                    if right_coord > 0 and left_coord > 0 and right_coord > left_coord:
                        outrow = '%s\t%s\t%s\t%s\t0\t%s\n' %(chromosome,str(left_coord),str(right_coord),ID,str(strand)) 
                        outputf.write(outrow)
                             
        inputf.close()
        outputf.close()
        elapsed = str(datetime.datetime.now() - self.starttime)
        print("There are %i warning(s)\n %s\n" %(len(self.warnings), '\n'.join(self.warnings)))
        print("Created file: %s" %self.oname)
        print('End date and time: %s\nElapsed time: %s\nDone!' %(str(datetime.datetime.now()),elapsed))
        sys.stdout.flush()
        


    
parser = OptionParser()
parser.add_option("-f", "--fname", action="store", type="string", dest="fname", metavar="<file>", help="Complete path to input file", default=-1)
parser.add_option("-o", "--oname", action="store", type="string", dest="oname", metavar="<file>", help="Complete path to output file", default=-1)
parser.add_option("-w", "--window", action="store", type="string", dest="window", metavar="<string>", help="[optional] Coordinates for window. E.g. promoter coordinates centred at +-1kb of the TSS --window=TSS:-1000:1000; Other examples --window=TES:0:2000; --window=GeneBody:-1000:1000. Can also be used to define windows relative to coordinates in any column, for instance, a 500bp window flanking the coordinates on column2 would be defined as  2:-500:500. Default is None", default=None)
parser.add_option("-t", "--sep", action="store", type="string", dest="sep", metavar="<string>", help="[optional] Separator of file (tab, comma, semicolon, space). Default is --sep=Tab", default='Tab')
parser.add_option("-c", "--columns", action="store", type="string", dest="columns", metavar="<string>", help="[optional] chr,start,end,strand columns of file separated by commas. Strand is optional but needed if --wstrand. Default is --columns=1,2,3", default='1,2,3')
parser.add_option("-i", "--cID", action="store", type="string", dest="cID", metavar="<string>", help="[optional] Column nr of the ID column (e.g.: --cID=1). Default is None", default=None)
parser.add_option("-s", "--wstrand", action="store_true", dest="wstrand", help='[optional] Define window coordinates based on strand. This takes into account the orientation of the feature to define "upstream" and "downstream" coordinates. Particularly important if deffining assymetrical windows', default=False)
parser.add_option("--onebased", "--onebased", action="store_true", dest="onebased", help='Specify this option if your input file has 1-based start positions. For instance, the GFF format from Ensemble uses 1-based coordinates for both the start and the end positions. While UCSC formats (such as Bed files, or other files downloaded from the Table Browser) are 0-based at the start and 1-based at the end positions. By default, this script uses start positions as if they were 0-based.', default=False)

(opt, args) = parser.parse_args(sys.argv)
if len(sys.argv) < 3: 
    parser.print_help()
    sys.exit(1)

for status in GenomicInterval(fname = opt.fname, oname = opt.oname, columns = opt.columns, sep = opt.sep,
                cID = opt.cID, window = opt.window, wstrand = opt.wstrand, onebased = opt.onebased).go(): print(status)

