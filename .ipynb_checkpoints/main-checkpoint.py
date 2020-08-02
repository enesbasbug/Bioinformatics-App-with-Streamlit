import streamlit as st
import neatbio.sequtils as utils
from Bio.Seq import Seq
from Bio import SeqIO
from collections import Counter

# data pkgs
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
import numpy as np
st.set_option('deprecation.showfileUploaderEncoding', False)

def main():

    st.title('Simple Bioinformatics App')
    menu = ["DNA Sequence","Dot Plot"]
    choice = st.sidebar.selectbox("Select Activity",menu)

    

    if choice == "DNA Sequence":
        st.subheader("DNA Sequence Analysis")
        seq_file = st.file_uploader("Upload FASTA File",type = ["fasta","fa"])
        #text_io = io.TextIOWrapper(seq_file)

        if seq_file is not None:
            dna_record = SeqIO.read(seq_file,"fasta")
            #st.write(dna_record)
            dna_seq = dna_record.seq
            desc = dna_record.description
            details = st.radio("Details", ("Description", "Sequence"))
            if details == "Description":
                 st.write(desc)
            elif details == "Sequence":
                 st.write(dna_seq)

            # Nucleotide Frequencies
            st.subheader("Nucleotide Frequencies")
            dna_freq = Counter(dna_seq)
            st.write(dna_freq)
            adenine_color = st.beta_color_picker("Adenine Color")
            guanine_color = st.beta_color_picker("Guanine Color")
            thymine_color = st.beta_color_picker("Thymine Color")
            cytosil_color = st.beta_color_picker("Cytosil Color")



            if st.button("Plot Freq"):
                barlist = plt.bar(dna_freq.keys(), dna_freq.values())
                barlist[0].set_color(adenine_color)
                barlist[1].set_color(thymine_color)
                barlist[2].set_color(guanine_color)
                barlist[3].set_color(cytosil_color)

                st.pyplot()


            st.subheader("DNA Composition")
            gc_score = utils.gc_content(str(dna_seq))
            at_score = utils.at_content(str(dna_seq))
            st.json({"GC Content ":gc_score, "AT Content ":at_score})

            # Nucleotide Count
            nt_count = st.text_input("Enter Nucleotide Here","Type Nucleotide Alphabet")
            st.write("Number of {} nucleotide is : {} ".format((nt_count),str(dna_seq).count(nt_count)))

            # Protein Synthesis
            p1 = dna_seq.translate()
            aa_freq = Counter(str(p1))

            st.subheader("Protein Synthesis")
            if st.checkbox("Transcription"):
                st.write(dna_seq.transcribe())
            elif st.checkbox("Translation"):
                st.write(dna_seq.translate())
            elif st.checkbox("Complement"):
                st.write(dna_seq.complement())
            elif st.checkbox("AA Frequency"):
                st.write(aa_freq)
            elif st.checkbox("AA Plot Frequency"):
                #aa_color = st.beta_color_picker("Amino Acid Color")
                #barlist = plt.bar(aa_freq.keys(), aa_freq.values(), color = aa_color)
                plt.bar(aa_freq.keys(), aa_freq.values())
                #barlist[0].set_color(aa_color)
                st.pyplot()
            elif st.checkbox("Full Amino Acid Name"):
                aa_name = str(p1).replace("*","")
                st.write(aa_name)
                st.write("--------------------------")
                st.write(utils.convert_1to3(aa_name))



    elif choice == "Dot Plot":
        st.subheader("Generate Dot Plot For Two Sequences")
        seq_file = st.file_uploader("Upload 1st FASTA File", type=["fasta", "fa"])
        seq_file2 = st.file_uploader("Upload 2nd FASTA File", type=["fasta", "fa"])

        # text_io = io.TextIOWrapper(seq_file)

        if seq_file and seq_file2 is not None:
            dna_record1 = SeqIO.read(seq_file, "fasta")
            dna_record2 = SeqIO.read(seq_file2, "fasta")

            # st.write(dna_record)
            dna_seq1 = dna_record1.seq
            dna_seq2 = dna_record2.seq

            desc1 = dna_record1.description
            desc2 = dna_record2.description

            details = st.radio("Details", ("Description", "Sequence"))
            if details == "Description":
                st.write(desc1)
                st.write("----------")
                st.write(desc2)
            elif details == "Sequence":
                st.write(dna_seq1)
                st.write("----------")
                st.write(dna_seq2)

            custom_limit = st.number_input("Select max number of Nucleotide ",10,200,25)
            if st.button("Dot Plot"):
                st.write("Comparing the first {} Nucleotide of Two Sequences ".format(custom_limit))
                dotplotx(dna_seq1[0:custom_limit], dna_seq2[0:custom_limit])
                st.pyplot()

    

def delta(x,y):
    return 0 if x == y else 1

def M(seq1,seq2,i,j,k):
    return sum(delta(x,y) for x,y in zip(seq1[i:i+k], seq2[j:j+k]))

def makeMatrix(seq1,seq2,k):
    n = len(seq1)
    m = len(seq2)
    return [[M(seq1,seq2,i,j,k) for j in range(m-k+1)] for i in range(n-k+1)]

def plotMatrix(M, t, seq1, seq2, nonblank = chr(0x25A0), blank = ' '):
    print(' |' + seq2)
    print('-'*(2 + len(seq2)))
    for label,row in zip(seq1, M):
        line = ''.join(nonblank if s < t else blank for s in row)
        print(label + '|' + line )

def dotplot(seq1, seq2, k=1, t=1):
    M = makeMatrix(seq1,seq2,k)
    plotMatrix(M,t,seq1,seq2)

    # Convert to Fxn
def dotplotx(seq1,seq2):
    plt.imshow(np.array(makeMatrix(seq1,seq2,1)))
    # on x axis list all sequnces of seq2
    xt = plt.xticks(np.arange(len(list(seq2))), list(seq2))
    # on y axis list all sequnces of seq1
    yt = plt.yticks(np.arange(len(list(seq1))), list(seq1))
    plt.show()



if __name__ == '__main__':
    main()