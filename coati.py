from Bio import Phylo
from Bio import SeqIO
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio import Entrez
from Bio import SeqIO
from Bio import GenBank 
import csv 
import re 
import matplotlib
import matplotlib.pyplot as plt
from Bio.Align.Applications import ClustalwCommandline
import os

##Primera función.- 
def fasta_downloader():
    '''
    La siguiente función con el nombre fasta_downloader nos permite ejecutar:
    Primero con el with open se convierte el archivo coati.txt en una lista, posterior a esto se hace uso del entrez de Biopythony con la
    finalidad de obtener la secuencia de cada identificador en el formato genbank y finalmente guardamos los datos en un documento .gb
    '''
    with open('data\coati.txt', 'r+') as coatlist:
        id_coati= coatlist.readlines()
        coati= []
        Entrez.email= "jozzfer14@gmail.com"
        for b in id_coati:
            handle = Entrez.efetch(db="nucleotide", id=b, rettype="gb", retmode="text") 
            seq_record= SeqIO.read(handle, "genbank")
            coati.append(seq_record)
            SeqIO.write(coati, "data\coati.gb", "genbank")
            
##Segunda función
def alignment():
    '''
    La siguiente función con el nombre alignment nos permite ejecutar:
     Primero transformar el archivo que se tiene de tipo genbank en formato fasta definiendo las variables que se presentan a continuación y posterior a esto se realizó un alineamiento de las secuencias empleando ClustalW con ayuda de la variable culstalw_exe y así poder acceder al programa ClustalW2, una vez realizado el alineamiento el resultado es guardado en los archivos coati.aln y coati.dnd.
     '''
    sequences= SeqIO.parse('data\coati.gb', 'genbank')
    SeqIO.write(sequences, "data\coati.fasta", "fasta")       
    clustalw_exe = r"C:\Archivos de programa (x86)\ClustalW2\clustalw2.exe"
    clustalw_cline = ClustalwCommandline(clustalw_exe, infile="data\coati.fasta")
    assert os.path.isfile(clustalw_exe),"Clustal_W executable is missing or not found"
    stdout, stderr = clustalw_cline()
    print(clustalw_cline)
    alignment = AlignIO.read(aln,"clustal")
    AlignIO.write(alignment, "data/coati.aln", "fasta")
    AlignIO.write(alignment, "data/coati.dnd", "phylip-relaxed")
    
##Tercera función
def tree():
    """
    La función con el nombre tree es capaz de realizar un cálculo de distancias haciendo uso del archivo aln para poder formar un árbol filogenético y su resultado poder guardarlo en un archivo pdf
    """
    with open("data\coati.aln", "r") as f:
        alignment  = AlignIO.read(f,"clustal")
        distance_calculation= DistanceCalculator("identity")
        dmatrix= distance_calculation.get_distance(alignment)
        tree_constructor=DistanceTreeConstructor(distance_calculation)
        tree=tree_constructor.build_tree(alignment)
        tree.rooted= True
        Phylo.draw(tree)
        plt.show()
        plt.savefig("data/coati_phylotree.pdf")
