

def match_peptide_to_protein(Protein_peptides,Reference_Proteome,cpus=1):
    # Here we determine all the protrein identities if not already provided in file.
    print("Assigning protein IDs based on search")
    Protein_peptides["Protein"]=""
    # sometimes data outputs contain a peptide sequence in brackets - the folowing will remove this
    Protein_peptides.Peptide=Protein_peptides.Peptide.str.replace(".*\]\.",'')
    Protein_peptides.Peptide=Protein_peptides.Peptide.str.replace("\.\[.*","")
    peptide_protein_identities = Determine_Peptide_Protein_Identities(Protein_peptides.Peptide.unique(),Reference_Proteome,cpus=cpus).determine_all_proteins()
    for key in peptide_protein_identities.keys():
        Protein_peptides.loc[Protein_peptides.Peptide == key,"Protein"]=peptide_protein_identities[key]
    return Protein_peptides

class Determine_Peptide_Protein_Identities:
    # This class takes in the peptide sequences
    # and determines the protein it may belong to
    def __init__(self,peptides,Reference_Proteome,cpus=1):
        self.peptides = peptides
        self.Reference_Proteome = Reference_Proteome
        self.protein_peptides={}
        self.cpus = cpus
        
    def append_results(self,result):
        # We append the results here
        self.protein_peptides.update(result)
        
    def retrieve_all_proteins(self,peptide,Reference_Proteome):
        # Here we find all the protein IDs that contain the peptide sequence
        Proteins_containing_peptide = Reference_Proteome[Reference_Proteome['FASTA'].str.contains(peptide)]["Uniprot_ID"]
        All_Proteins = ";".join(Proteins_containing_peptide)
        return {'peptide':peptide,'All_Proteins':All_Proteins}
    
    def find_matches(self,to_process):
        all_peptides={}
        for peptide in to_process:
            results = self.retrieve_all_proteins(peptide,self.Reference_Proteome)
            all_peptides[results['peptide']]=results['All_Proteins']
        return all_peptides
        
    def determine_all_proteins(self):
        import multiprocessing as mp  
        pool = mp.Pool(self.cpus)
        count= 0
        # Here we split the peptides in batches
        to_process=[]
        for peptide in self.peptides:
            count+=1
            if (count % 1000 == 0):
                if self.cpus>1:
                    pool.apply_async(self.find_matches,args=([to_process]),callback=self.append_results)
                else:
                    result = self.find_matches(to_process)
                    self.append_results(result)
                to_process=[]
            else:
                # Here we append to the list to process
                to_process.append(peptide)
        # Here we process the remaining ones that dint get triggered by exacution 
        if self.cpus>1:
            pool.apply_async(self.find_matches,args=([to_process]),callback=self.append_results)
        else:
            result = self.find_matches(to_process)
            self.append_results(result)
        pool.close()
        pool.join() 
        return self.protein_peptides
