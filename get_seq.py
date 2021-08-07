# UPDATE 01/08/2020: 
######## length leveled for b, l, and r.
from Bio import SeqIO
import argparse
from argparse import RawTextHelpFormatter
from itertools import product 
parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument("-f", "--fasta", help = "input fasta file", required=True)
parser.add_argument("-p", "--positions", help = "a file that contains list of positions to extract flanking sequences", required=True)
parser.add_argument("-m", "--mode", help = "which flanking, r=right, l=left, b=both", type=str, choices=["r", "l", "b"], required=True)
parser.add_argument("-l", "--length", help = "length of the flanking sequence", type = int, required=True)
parser.add_argument("-v", "--vcf", help ="input vcf file (optional)", required = False)
args = parser.parse_args()


seqs = list(SeqIO.parse(args.fasta, "fasta"))
positions = open (args.positions, "r")
key=[]
key1=[]
value=[]
value1=[]

def get_chr(c):
    if c=="X":
        return 23
    elif c=="Y":
        return 24
    elif c=="M":
        return 25
    else:
        return int(c)


def replacer(s, newstring, index, nofail=False):
    # raise an error if index is outside of the string
    if not nofail and index not in range(len(s)):
        raise ValueError("index outside given string")

    # if not erroring, but the index is still not in the correct range..
    if index < 0:  # add it to the beginning
        return newstring + s
    if index > len(s):  # add it to the end
        return s + newstring

    # insert the new string between "slices" of the original
    return s[:index] + newstring + s[index + 1:]

def find_possible_variant_sequences(refSequence, variants):
  """ Find all possible variant sequences
  
  Arguments:
      refSequence {str} -- reference sequence
      variants {List of array} -- List of array with location (zero-based) and possible bases, for examples [[3, 'C', 'T'], [10, 'C', 'G', 'T']]
  """
  resSeq = list(refSequence)
  variantsBases = [v[1:] for v in variants]
  res = list(product(*variantsBases)) 

  result = []
  for resValues in res:
    for i in range(len(variants)):
      variant = variants[i]
      resSeq[variant[0]] = resValues[i]
    result.append("".join(resSeq))
        
  return(result)


output = open("human_"+args.mode+"_flanking_"+str(args.length)+".txt", 'w+')
# no VCF file
if args.vcf== None:
    if args.mode == "r":
    #output = open("human_"+args.mode+"_flanking_"+str(args.length)+".txt", 'w')
        for line in positions:
            words=line.split("\t")
            chr=get_chr(words[0])-1
            p=int(words[1])
            s=str(seqs[chr].seq[p:p+args.length])
            output.write (line.strip("\n").strip("\r")+"\t"+s+"\n")
        output.close()
    
    if args.mode == "l":
    #output = open("human_"+args.mode+"_flanking_"+str(args.length)+".txt", 'w')
        for line in positions:
            words=line.split("\t")
            chr=get_chr(words[0])-1
            p=int(words[1])
            s=str(seqs[chr].seq[p-args.length-1:p-1])
            output.write (line.strip("\n").strip("\r")+"\t"+s+"\n")
        output.close()


    if args.mode == "b":
    #output = open("human_"+args.mode+"_flanking_"+str(args.length)+".txt", 'w')
        for line in positions:
            words=line.split("\t")
            chr=get_chr(words[0])-1
            p=int(words[1])
            offset=args.length/2
            s=str(seqs[chr].seq[p-offset-1:p+offset]) # Yan used to write "[p-offset-1:p+offset]"
            output.write (line.strip("\n").strip("\r")+"\t"+s+"\n")
        output.close()

#has VCF file
else:  
    #output = open("human_"+args.mode+"_flanking_"+str(args.length)+".txt", 'w')
    #extract all the positions information out
    if args.mode == "r":
        with open(args.vcf) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                else:
                    line = line.strip()
                    line= line.replace('"','')
                    line= line.replace('rs','')
                    line= line.split("\t") 
                    if ('A' in line[3] or 'C' in line[3] or 'G' in line[3] or 'T' in line[3] or '.' in line[3]) and ('A' in line[4] or 'C' in line[4] or 'G' in line[4] or 'T' in line[4]):                     
                # Add SNP positions
                        key.append(line[1])
                            #key.append(line[2])
                            #key.append(line[3])
                # Add SNP informations       
                        value.append([line[1],line[3],line[4]])
                            #value.append([line[2],line[4],line[5]])
                            #value.append([line[3],line[4],line[5]]) 
                    else:
                        if "." not in line[5]:
                            key.append(line[1])
                            value.append([line[1],line[4],line[5]])

#separate positions that has SNP
        for line in positions:
            a=1
            line = line.strip()
            words=line.split("\t")
            chr=get_chr(words[0])-1
            p1 = words[1]
            p=int(words[1])
            # If the postions no in vcf files
            if p1 not in key:
                a=1
                position1=[]
                s=str(seqs[chr].seq[p:p+args.length])
                output.write(line.strip("\n").strip("\r")+"\t"+s+"\n")
                while (a<=args.length):                   
                    if str(p+a) in key:                       
                        index2= key.index(str(p+a))
                        if value[index2][1] =='.':
                            value[index2][1] = s[0+int(value[index2][0])-p].upper()
                        if type(value[index2][2]) ==str:
                            value[index2][2]=value[index2][2].split(',')
                        if len(value[index2][2])==1:
                            position1.append([(int(value[index2][0])-p),value[index2][1],value[index2][2][0]])
                        else:
                            position1.append([(int(value[index2][0])-p),value[index2][1],value[index2][2][0],value[index2][2][-1]])
                    a=a+1
                #print position1
                result = find_possible_variant_sequences(s,position1)
                for i in result:
                    i = i[0:args.length]
                    output.write(line.strip("\n").strip("\r")+"\t"+i+"\n")

            # If the positions in vcf files (insertion and replacement)
            elif p1 in key: 
                a=1             
                position=[]
                index1 = key.index(p1)
            	s = str(seqs[chr].seq[p:p+args.length])
            	# Get the index of the position in the key list
                if value[index1][1] =='.':
                    value[index1][1] = s[0].upper()
                if type(value[index1][2]) ==str:
                    value[index1][2]=value[index1][2].split(',')
                if len(value[index1][2])==1:
                    position.append([(int(value[index1][0])-p),value[index1][1],value[index1][2][0]])
                else:
                    position.append([(int(value[index1][0])-p),value[index1][1],value[index1][2][0],value[index1][2][-1]])            
                while (a<=args.length):                   
                    if str(p+a) in key:                       
                        index2= key.index(str(p+a))
                        if value[index2][1] =='.':
                            value[index2][1] = s[0+int(value[index2][0])-p].upper()
                        if type(value[index2][2]) ==str:
                            value[index2][2]=value[index2][2].split(',')
                        if len(value[index2][2])==1:
                            position.append([(int(value[index2][0])-p),value[index2][1],value[index2][2][0]])
                        else:
                            position.append([(int(value[index2][0])-p),value[index2][1],value[index2][2][0],value[index2][2][-1]])
                    a=a+1
                result = find_possible_variant_sequences(s,position)
                for i in result:
                    i = i[0:args.length]
                    output.write(line.strip("\n").strip("\r")+"\t"+i+"\n")
            # If the positions in vcf files (deletion)

        output.close()

    if args.mode == "l":
        with open(args.vcf) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                else:
                    line = line.strip()
                    line= line.replace('"','')
                    line= line.replace('rs','')
                    line= line.split("\t") 
                    if ('A' in line[3] or 'C' in line[3] or 'G' in line[3] or 'T' in line[3] or '.' in line[3]) and ('A' in line[4] or 'C' in line[4] or 'G' in line[4] or 'T' in line[4]) :                                                    
                # Add SNP positions
                        key.append(line[1])
                            #key.append(line[2])
                            #key.append(line[3])
                # Add SNP informations       
                        value.append([line[1],line[3],line[4]])
                            #value.append([line[2],line[4],line[5]])
                            #value.append([line[3],line[4],line[5]]) 
                    else:
                        if "." not in line[5]:
                            key.append(line[1])
                            value.append([line[1],line[4],line[5]])
                   
#separate positions that has SNP
        for line in positions:
            line = line.strip()
            words=line.split("\t")
            chr=get_chr(words[0])-1
            p1 = words[1]
            p=int(words[1])
            # If the postions no in vcf files
            if p1 not in key:
                a=1
                position1=[]
                #s=str(seqs[chr].seq[p:p+args.length])
                output.write(line.strip("\n").strip("\r")+"\t"+s+"\n")
                while (a<=args.length):                   
                    if str(p-a) in key:                       
                        index2= key.index(str(p+a))
                        if value[index2][1] =='.':
                            value[index2][1] = s[0+int(value[index2][0])-p].upper()
                        if type(value[index2][2]) ==str:
                            value[index2][2]=value[index2][2].split(',')
                        if len(value[index2][2])==1:
                            position1.append([(int(value[index2][0])-p),value[index2][1],value[index2][2][0]])
                        else:
                            position1.append([(int(value[index2][0])-p),value[index2][1],value[index2][2][0],value[index2][2][-1]])
                    a=a+1
                #print position1
                result = find_possible_variant_sequences(s,position1)
                for i in result:
                    i = i[-args.length:]
                    output.write(line.strip("\n").strip("\r")+"\t"+i+"\n")
            # If the positions in vcf files (insertion and replacement)
            elif p1 in key: 
                a=1              
                position=[]
                index1 = key.index(p1)
                s=str(seqs[chr].seq[p-args.length-1:p-1])
                # Get the index of the position in the key list
                if value[index1][1] =='.':
                    value[index1][1] = s[-1].upper()
                if type(value[index1][2]) ==str:
                    value[index1][2]=value[index1][2].split(',')
                if len(value[index1][2])==2:
                    position.append([(int(value[index1][0])-p-1),value[index1][1],value[index1][2][0]])
                else:
                    position.append([(int(value[index1][0])-p-1),value[index1][1],value[index1][2][0],value[index1][2][1]])  

                while (a<=args.length):                   
                    if str(p-a) in key:                       
                        index2= key.index(str(p-a))
                        if value[index2][1] =='.':
                            value[index2][1] = s[0+int(value[index2][0])-p].upper()
                        if type(value[index2][2]) ==str:
                            value[index2][2]=value[index2][2].split(',')
                        if len(value[index2][2])==2:
                            position.append([(int(value[index2][0])-p-1),value[index2][1],value[index2][2][0]])
                        else:
                            position.append([(int(value[index2][0])-p-1),value[index2][1],value[index2][2][0],value[index2][2][1]])
                    a=a+1
        
                result = find_possible_variant_sequences(s,position)

                for i in result:
                    i = i[-args.length:]
                    output.write(line.strip("\n").strip("\r")+"\t"+i+"\n")
            # If the positions in vcf files (deletion)
        output.close()

    if args.mode == "b":
        with open(args.vcf) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                else:
                    line = line.strip()
                    line= line.replace('"','')
                    line= line.replace('rs','')
                    line= line.split("\t") 
                    if ('A' in line[3] or 'C' in line[3] or 'G' in line[3] or 'T' in line[3] or '.' in line[3]) and ('A' in line[4] or 'C' in line[4] or 'G' in line[4] or 'T' in line[4]) :                  
                                                            
                # Add SNP positions
                            key.append(line[1])
                            #key.append(line[2])
                            #key.append(line[3])
                # Add SNP informations       
                            value.append([line[1],line[3],line[4]])
                            #value.append([line[2],line[4],line[5]])
                            #value.append([line[3],line[4],line[5]]) 
                    else:
                        if "." not in line[5]:
                            key.append(line[1])
                            value.append([line[1],line[4],line[5]])

        for line in positions:
            line = line.strip()
            words=line.split("\t")
            chr=get_chr(words[0])-1
            offset=args.length/2
            p1 = words[1]
            p=int(words[1])
            # If the postions no in vcf files
            if p1 not in key:
                a=1
                b=1
                #print "postions ", p
                position1=[]
                s=str(seqs[chr].seq[p-offset-1:p+offset])
                #output.write(line.strip("\n").strip("\r")+"\t"+s+"\n")                
                while (a<=offset):                   
                    if str(p-a) in key:                       
                        index2= key.index(str(p-a))
                        if value[index2][1] =='.':
                            value[index2][1] = s[0+int(value[index2][0])-p+offset].upper()
                        if type(value[index2][2]) ==str:
                            value[index2][2]=value[index2][2].split(',')
                        if len(value[index2][2])==2 or len(value[index2][2])==1:
                            position1.append([(int(value[index2][0])-p+offset),value[index2][1],value[index2][2][0]])
                        else:
                            position1.append([(int(value[index2][0])-p+offset),value[index2][1],value[index2][2][0],value[index2][2][1]])
                    a=a+1
                

                while (b<=offset):
                    if str(p+b) in key:
                        #print p+b
                        index3= key.index(str(p+b))
                        if value[index3][1] =='.':
                            value[index3][1] = s[0+int(value[index3][0])-p+offset].upper()
                        if type(value[index3][2]) ==str:
                            value[index3][2]=value[index3][2].split(',')
                        if len(value[index3][2])==2 or len(value[index3][2])==1:
                            position1.append([(int(value[index3][0])-p+offset),value[index3][1],value[index3][2][0]])
                        else:
                            position1.append([(int(value[index3][0])-p+offset),value[index3][1],value[index3][2][0],value[index3][2][1]])
                    b=b+1
                result = find_possible_variant_sequences(s,position1)
                # Pay attention to the length of the sequence.
                for i in result:
                    i = i[offset-args.length/2:offset]+i[offset]+i[offset+1:offset+args.length/2+1]
                    output.write(line.strip("\n").strip("\r")+"\t"+i+"\n")

            # If the positions in vcf files (insertion and replacement)
            elif p1 in key:  
                a=1
                b=1           
                position=[]
                index1 = key.index(p1)
                s=str(seqs[chr].seq[p-offset-1:p+offset])
                # Get the index of the position in the key list
                if value[index1][1] =='.':
                    value[index1][1] = s[offset].upper()
                if type(value[index1][2]) ==str:
                    value[index1][2]=value[index1][2].split(',')
                if len(value[index1][2])==2 or len(value[index1][2])==1:
                    position.append([(int(value[index1][0])-p+offset),value[index1][1],value[index1][2][0]])
                else:
                    position.append([(int(value[index1][0])-p+offset),value[index1][1],value[index1][2][0],value[index1][2][1]])  

                while (a<=offset):                   
                    if str(p-a) in key:                       
                        index2= key.index(str(p-a))
                        if value[index2][1] =='.':
                            value[index2][1] = s[0+int(value[index2][0])-p+offset].upper()
                        if type(value[index2][2]) ==str:
                            value[index2][2]=value[index2][2].split(',')
                        if len(value[index2][2])==2 or len(value[index2][2])==1:
                            position.append([(int(value[index2][0])-p+offset),value[index2][1],value[index2][2][0]])
                        else:
                            position.append([(int(value[index2][0])-p+offset),value[index2][1],value[index2][2][0],value[index2][2][1]])
                    a=a+1
                while (b<=offset):
                    if str(p+b) in key:
                        index3= key.index(str(p+b))
                        if value[index3][1] =='.':
                            value[index3][1] = s[0+int(value[index3][0])-p+offset].upper()
                        if type(value[index3][2]) ==str:
                            value[index3][2]=value[index3][2].split(',')
                        if len(value[index3][2])==2 or len(value[index3][2])==1:
                            position.append([(int(value[index3][0])-p+offset),value[index3][1],value[index3][2][0]])
                        else:
                            position.append([(int(value[index3][0])-p+offset),value[index3][1],value[index3][2][0],value[index3][2][1]])
                    b=b+1
                result = find_possible_variant_sequences(s,position)
                # Pay attention to the length of the sequence.
                for i in result:
                    i = i[offset-args.length/2:offset]+i[offset]+i[offset+1:offset+args.length/2+1]
                    output.write(line.strip("\n").strip("\r")+"\t"+i+"\n")
        output.close()    




                    


         
    