from Bio import Entrez
import os



# Configuração do email (obrigatório para usar o Entrez)
Entrez.email =  "lucasdias26@hotmail.com.br"

# String de busca
query = '"Severe acute respiratory syndrome coronavirus 2"[Organism] OR sars-cov-2[All Fields]) AND complete[All Fields] AND genome[All Fields] NOT mRNA[All Fields] NOT RNA[All Fields'

# 1. Realiza a busca no banco de dados "nucleotide"
handle = Entrez.esearch(db="nucleotide", term=query, retmax=20)  # retmax define o número máximo de resultados, o padão é 20
search_results = Entrez.read(handle)
handle.close()

# Obtém os IDs dos resultados
id_list = search_results["IdList"]
print(f"Total de genomas encontrados: {len(id_list)}")

# 2. Baixa os genomas completos no formato FASTA
if id_list:
    handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="fasta", retmode="text")
    genomas = handle.read()
    handle.close()

    # Salva os genomas em um arquivo FASTA
    here = os.path.dirname(os.path.abspath(__file__))
    destinePath = os.path.join(here, "genomas_sars_cov_2.fasta")
    
    with open(destinePath, "w") as file:
        file.write(genomas)
    print("Genomas baixados e salvos em 'genomas_sars_cov_2.fasta'.")
else:
    print("Nenhum genoma encontrado.")