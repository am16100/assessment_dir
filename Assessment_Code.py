

from Bio.Seq import Seq 
from Bio import SeqIO #used biopython as the desired functions are already in the package and so requires less lines of code#

#first take sequence from gen bank that you want to design primers for - have to ensure no spaces
coding_dna = Seq('ATGGCTATGAAACAAAGCAAAACGGTGTACATGTTTCCCGGGCAAGGCTCGCAGTTTCGCGGAATGGGTGAGGGCTTGTTCGAACGCTTTGCCGAACTGACGGCCTGTGCCGATCGCGTGCTGGGCTACTCGATCCGTGAGCTGTGCGAAAACGATCCGCGCAACGAACTTGGCCAGACCCGCTTCACCCAACCTGCCCTGTATGTGGTCAATGTACTGTCGTACCTGGCTCAGGCTGGCGATGCCGCGCCCCCCGACTATGTGCTCGGGCACAGCCTCGGCGAGTTCTGCGCTTTGTTTGCCGCCGGCGCCTATGACTTCGAAACCGGATTGCGCCTGGTCAAGCGGCGCGGCGAACTCATGTCCGAGGCCACCGGGGGCGGCATGTCGGCTGTGCTGAATCTGGATCTTGCGACGATCAAACAGGTGTTGCGCCAGGCTGGGAGTACACAGCTGGATTTCGCCAACTTCAATGCCCCGCAACAGACCGTATTGGCCGGACCGCTCGACGCACTGGAAAGCGTGCGTACGCAGATTGAAGACGCCTCGGGTATTTGTGTCGCGTTGAATGTCAGCGCGCCTTTCCACTCCCGCTACATGCGTGGCGCCCAAGAAGCGTTTGCGGCTGAGCTGGCCCGTGTGACCTTCAAGCCCCTGACGCTGCCGGTGATCGCCAATGTCGACGCACGCCCTTATGAGCAGGAAGCTATCGCCAGCCAGTTGGCGCGGCAAATGACGTCATCGGTGCAGTGGGTCGAAAGCATCGAGTATTTACTGCAGGCCGGCATCACCCAGTTCAAGGAGATCGGCCCCGGCAATGTATTGACCAATCTGCAAGCGAAAATCGAAAAAAACCGTTCGCCAGCCCAGCCAGTCGTGGCATCTGCCGTGGTTGCACCGCTGCCTTTGGCTGCGCGGCATGACGACGGCGTGCCGGGCCTGGCGCTGACAGCCGAGGGGTTGGGGAGCGCC')

    
#design primers, requires the adapter sequence to be added to the first 20 bases and the last 20
    
popinf_primer_for_adapter = Seq('AAGTTCTGTTTCAGGGCCCG')

popinf_primer_rev_adapter = Seq('ATGGTCTAGAAAGCTTTA')

for_primer = coding_dna[0:21] #defines region of the DNA sequence needed to create the primer

popF_for_primer = popinf_primer_for_adapter + for_primer #adds the adapter sequence to the short section of the DNA sequence
print(f"Forward primer sequence:{popF_for_primer}") # output the primer sequence to be produced - added some text to clarify what each output is when the code is run in the terminal
    
rev_primer = (coding_dna[-20:])[::-1] # managed to get the last 20 but was the wrong direction so added [::-1] to reverse the sequence
popF_rev_primer = popinf_primer_rev_adapter + rev_primer
print(f"Reverse primer sequence:{popF_rev_primer}")
    
#QC check the primer design#   
from Bio.SeqUtils import GC
GC_for = GC(popF_for_primer)
GC_rev = GC(popF_rev_primer) # gives the GC content but this in itself is not useful so added a if else loop
    
GC_diff = abs(GC_for-GC_rev) #wanted it to always be positive so added absolute function
if GC_diff <=20:
    print('GC difference: pass')
    
else:
    print('GC difference: fail')
    
messenger_rna = coding_dna.transcribe() #initally tried to use a dictionary but this did not work therefore used biopython
    
my_protein = messenger_rna.translate(table="Bacterial") # added in the table argument as this DNA sequence was taken from bacteria
print(f"Protein sequence:{my_protein}")

import MW                                #wrote the function but then decided to create a module as this code could be reused in various other contexts
Molecular_weight = MW.find_MW(my_protein)

print(f"MW = {Molecular_weight}")