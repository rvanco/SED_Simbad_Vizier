"""
Created on Wed Nov 16 10:01:00 2022

@author: Renaud Vancoellie
"""

import Def as D

path = "result/"

variation = int(input("0 - Use and Download all the question \n1 - Choose the Question \n[Which option would you like ?] "))

    
if variation == 0 :
    print("\nSimbad target resolver :")
    tbl_coord, tbl_flux = D.simbad_target_name_resolver(path)
    print("\n\n\n\n\n\nVizier Cone search :")
    tbl_, tbl_coord, nb_filter, lambda_list, radius = D.vizier_cone_search(path)
    print("\n\n\n\n\n\nSED (Using Vizier and Simbad) :")
    D.SED(path, 2)
    
else :
    finish = "n"
    
    while finish != "y" :
    
        Choose = int(input("0 - Simbad target resolver \n1 - Vizier cone search \n2 - SED (using Vizier adn Simbad) \n[What function do you want to use ?] "))
        
        if Choose == 0 :
            print("\n\nSimbad target resolver :")
            tbl_coord, tbl_flux = D.simbad_target_name_resolver(path)
                
        elif Choose == 1 :
            print("\n\nVizier Cone search :")
            tbl_, tbl_coord, nb_filter, lambda_list, dlambda_list, radius = D.vizier_cone_search(path)
  
        else :
            print("SED (Using Vizier and Simbad) :")
            D.SED(path, 2)
            

        finish = input("Have you finish ? (y/n) ")
