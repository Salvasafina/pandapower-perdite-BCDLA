import numpy as np
import pandapower as pp
from funzioni_rete import *

'''
Si prova il metodo BCDLA:
- A partire dal file .mat che contiene l'informazione della rete si crea la rete in pandapower, questa rete ha 
un unico carico o generatore per nodo.
- Si fa un flusso di potenza con il metodo BWFS
- Si ricalcolano gli elementi usati nel BWSf necessari per applicare BCDLA: matrice di incidenze, matrice gamma,
 correnti ai rami, correnti ai nodi. Si devono calcolare perché il metodo BWFS di pandapower non li rende disponibili.
- Si applica il metodo BCDLA e si ottengono anche le perdite totali.
'''

# --- Creazione dei DF ausiliari e costruzione della rete in  Pandapower --- # 
# Chiamata alla funzione che crea la rete per Pandapower a partire del file .mat
#       Si crea un DF dei rami, un DF dei carichi, e un DF delle grandeze base per agevolare la costruzione della rete in pandapower
#       Crea una rete di pandapower vuota e aggiunge: nodi, rete esterna, linee, e carichi secondo i DF creati prima.
rete = crea_rete_pp('Rural_network.mat')
rete_df = dati_df_rete('Rural_network.mat')
print(rete)
# -------------------------------------------------------------------------- #


# ------- Soluzione della rete ------- # 
pp.runpp(rete, algorithm = 'bfsw')
# ------------------------------------ #

# ---- Ottenzione dati necessari per aplicare BCDLA ---- #
mat_incid = mat_incid(rete)
mat_gamma =np.linalg.inv(mat_incid)


i_long = correnti_rami(rete.line,rete.res_line, rete_df)
i_nodi = correnti_nodi(mat_gamma,i_long)
imp_real=resistenza_long(rete.line, rete_df)
# ------------------------------------------------------ #

# ---- Calcolo delle perdite con BCDLA ---- #
perdite =(i_nodi.conjugate()*(np.dot(mat_gamma,(imp_real * i_long)))).real
perd_tot = sum(perdite)

print(f'Perdite per nodo:\n {perdite}')
print(f'Perdite totali: \n{perd_tot}')
# ----------------------------------------- #