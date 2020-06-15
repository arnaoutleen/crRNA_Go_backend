from Bio.Seq import Seq

#use biopython library to 
string="ACTGGGGGGGGGGGGGGGGGGGGG"
input_seq=Seq(string)
# input is 3' to 5'. Then, just take its reverse complement using biopython command
new_seq=input_seq.reverse_complement()

#if the input_seq is a DNA then do this:
#new_seq=input_seq.transcribe().reverse_complement()

#print out the answer
result="3'-ACAACAUUAUCGGGGGUUUUGACCUGGAAGGUGUUG"+new_seq+"-5'"

print(result)