#added function into a module 

amino_acids = {'A':71, 'R':156, 'N':114, 'D':115, 'C':103, 'E':129, 'Q':128, 'G':57, 'H':137, 'I':113, 'L':113, 'K':128, 'M':131, 'F':147, 'P':97, 'S':87, 'T':101, 'W':186, 'Y':163,'V':99}    #created dictionary using the MW of the most comon amino acids

def find_MW(my_protein):
    aa = []
    for letter in my_protein:
        MW = amino_acids[letter]
        aa.append(MW)
        molecular_weight = sum(aa)   #the list aa only contained a list of MWs, had to sum the list
    return molecular_weight          #had to change the notation, originally haded the sum function after the return line but this would not work when the function was placed into a module.