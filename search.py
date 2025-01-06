from Bio import Entrez
import os



# Configuração do email (obrigatório para usar o Entrez)
Entrez.email = "lucasdias26@hotmail.com.br"  # Substitua pelo seu email

# String de busca
query = '"Severe acute respiratory syndrome coronavirus 2"[Organism] AND "complete genome"[All Fields] AND "human"[Host]'

# 1. Realiza a busca no banco de dados "nucleotide"
handle = Entrez.esearch(db="nucleotide", term=query, retmax=10)  # retmax define o número máximo de resultados
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