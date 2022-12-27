"""
Created on Wed Nov 16 09:57:06 2022

@author: rvancoellie
"""

########################################################################################################################################################################################################
            # IMPORT :
########################################################################################################################################################################################################

import urllib.request as http
from astropy.table import Table
from astropy import units as u
import numpy as np
import matplotlib.pyplot as plt
from cython_code import c_to_jsky as to_jsky
from cython_code import c_search_vega_filter as search_vega_filter

# Jansky Unit definition :
unit_jsky = u.def_unit('jsky')
nan = float("NaN")

########################################################################################################################################################################################################
            # SIMBAD TARGET RESOLVER :
########################################################################################################################################################################################################

def simbad_target_name_resolver(path, name="") :
    if name == "" :
        TARGET = input('[Object name you want to target] ')
    else :
        TARGET = name

    url = f"https://simbad.u-strasbg.fr/simbad/sim-id?output.format=ascii&Ident={TARGET}"
        
    nan = float("NaN")
    system = "Johnson"
    
    ar_filter_name = []
    ar_flux = []
    ar_flux_error = []
    ar_lambda = []
    ar_lambda_error = []
    
    tbl_coord = Table([[nan]*u.deg, [nan]*u.deg], names=("Ra ","Dec "), masked=True)

    with http.urlopen(url) as fd :
        for line in fd:
            if "Coordinates(ICRS,ep=J2000,eq=2000)" in line.decode('UTF-8'):
                r_line = line.decode('utf-8').strip().split()
                Ra_deg = ((float(r_line[1])  + float(r_line[2])/60 + float(r_line[3])/3600)/24) * 360
                Dec_deg = float(r_line[4])  + float(r_line[5])/60 + float(r_line[6])/3600
                tbl_coord = Table([[Ra_deg]*u.deg, [Dec_deg]*u.deg], names=("Ra ","Dec "))
                
            elif 'Flux' in line.decode('utf-8') :
                r_line = line.decode('utf-8').strip().split()
                
                if len(r_line) > 4 :
                    filter_name = r_line[1]
                    lambda_, dlambda, Fmag = search_vega_filter(system, filter_name)
                    
                    if Fmag == 0.0 :
                        flux = nan
                    else :
                        mag = float(r_line[3])
                        flux = to_jsky(Fmag, mag)

                    if lambda_ == 0.0 :
                        lambda_, dlambda = nan, nan
                    
                        
                    if len(r_line[4]) > 3 :
                        mag_neg = mag - float(r_line[4][1:-1])
                        flux_neg = to_jsky(Fmag, mag_neg)
                        mag_pos = mag + float(r_line[4][1:-1])
                        flux_pos = to_jsky(Fmag, mag_pos)
                        
                        error_neg = round(flux - flux_neg, 4)
                        error_pos = round(flux - flux_pos, 4)
                        
                        flux_error = f"+{error_pos}, {error_neg}"
                        
                    else :
                        print(f"WARNING : There is no flux error in Simbad data for the filter {filter_name}")
                        flux_error = f"+{nan}, -{nan}"
                    
                    ar_filter_name.append(filter_name)
                    ar_flux.append(flux)
                    ar_flux_error.append(flux_error)
                    ar_lambda.append(lambda_)
                    ar_lambda_error.append(dlambda)

    if tbl_coord["Ra "][0]*0 != 0 :
        print("[ERROR : Name selected is not valid.]")
    
                    
    tbl_flux = Table([ar_filter_name, ar_flux*unit_jsky, ar_flux_error, ar_lambda*u.um ,ar_lambda_error*u.um], 
                  names=("filter","flux","flux error (jsky)","Lambda","Lambda error"))
    
    tbl_coord.write(f'{path}tbl_coord_simbad.ecsv', overwrite=True)
    tbl_flux.write(f'{path}tbl_flux_simbad.ecsv', overwrite=True)
    
    return tbl_coord, tbl_flux




########################################################################################################################################################################################################
            # VIZIER CONE SEARCH :
########################################################################################################################################################################################################


def all_filter(i, lambda_list, dlambda_list, filter_list, column_list, system_list, table, tbl_):

    lambda_list.append([])
    dlambda_list.append([])
    if len(system_list) < 1 :
        system = system_list[i]
    else :
        system = system_list[0]
    filter_ = filter_list[i]
    column = column_list[i]       
    print(f"[The column {column} use the system {system} with the filter {filter_}.]")
    lambda_, dlambda, Fmag = search_vega_filter(system, filter_)
    flux_jsky_list = []
    err_flux_jsky_list = []

    try :
        column_err = "e_"+column
        try_table = table[column_err]
    except :
         column_err = "e"+column

    for j in range(len(table[column])) :
        table[column].fill_value = nan
        table[column_err].fill_value = nan
                
        mag_err = table[column_err].filled()[j]
        mag = table[column].filled()[j]
                
        mag_neg = mag - mag_err
        mag_pos = mag + mag_err
        flux = to_jsky(Fmag,mag)
        error_neg = round(flux - to_jsky(Fmag, mag_neg), 4)
        error_pos = round(flux - to_jsky(Fmag, mag_pos), 4)
        flux_error = f"+{error_pos}, {error_neg}"
                
        flux_jsky_list.append(flux)
        err_flux_jsky_list.append(flux_error)
                
        lambda_list[i].append(lambda_)
        dlambda_list[i].append(dlambda)
        
    tbl_.add_column(flux_jsky_list*unit_jsky, name=f"{column}")
    tbl_.add_column(err_flux_jsky_list, name=f"{column_err} (jsky)")
            
    for j in range(len(flux_jsky_list)) :
        if np.isnan(flux_jsky_list[j]) == True :
            print(f"[WARNING : There are some empty data in the filter {filter_}.]\n")
            break
                
    return tbl_

def while_filter(lambda_list, dlambda_list, filter_list, column_list, system_list, table, tbl_, filter_choose, iteration) :

    lambda_list.append([])
    dlambda_list.append([])
    column_nb = int(input(f"[Filter number {round(iteration,1)}] "))
    for i in range(len(filter_choose)) :
        if column_nb == filter_choose[i] :
            print(f"[ERROR : Do not select multiple time the same filter.]\nThe filter already selected are {filter_choose}")
            column_nb = int(input(f"[Filter number {round(iteration,1)}] "))
    try :
        system = system_list[column_nb]
        filter_ = filter_list[column_nb]
        column = column_list[column_nb]
    except :
        print(f"[ERROR : Choose a valid number (between 0 and {len(filter_list)})]")
        column_nb = int(input(f"[Filter number {round(iteration,1)}] "))
    filter_choose.append(column_nb)        
    print(f"\n[The column {column} use the system {system} with the filter {filter_}.]")
    lambda_, dlambda, Fmag = search_vega_filter(system, filter_)
    flux_jsky_list = []
    err_flux_jsky_list = []

    try :
        column_err = "e_"+column
        try_table = table[column_err]
    except :
        column_err = "e"+column

    for i in range(len(table[column])) :
        table[column].fill_value = nan
        table[column_err].fill_value = nan
             
        mag_err = table[column_err].filled()[i]
        mag = table[column].filled()[i]
                
        mag_neg = mag - mag_err
        mag_pos = mag + mag_err
        flux = to_jsky(Fmag,mag)
        error_neg = round(flux - to_jsky(Fmag, mag_neg), 4)
        error_pos = round(flux - to_jsky(Fmag, mag_pos), 4)
        flux_error = f"+{error_pos}, {error_neg}"
                
        flux_jsky_list.append(flux)
        err_flux_jsky_list.append(flux_error)
               
        lambda_list[iteration-1].append(lambda_)
        dlambda_list[iteration-1].append(dlambda)
                
    tbl_.add_column(flux_jsky_list*unit_jsky, name=f"{column}")
    tbl_.add_column(err_flux_jsky_list, name=f"{column_err} (jsky)")
            
    for j in range(len(flux_jsky_list)) :
        if np.isnan(flux_jsky_list[j]) == True :
            print(f"[WARNING : There are some empty data in the filter {filter_}.]")
            break
                
    return tbl_, filter_choose





def vizier_cone_search(path, filtre_SED="",catalogue="", center_name="", radius_orig="") :
    if catalogue == "" :
        name = input('0 - III/284/allstars \n1 - II/340/xmmom2_1 \n2 - II/262 \nFor other, write the name \n[number of the catalogue] ')
    else :
        name = str(catalogue)
    if center_name == "" :
        center = input("[Conesearch's center name (e.g. M77, HD1, NGC2264)] ")
    else :
        center = center_name
    if radius_orig == "" :
        radius = input("[Conesearch's radius (arcmin)] ")
    else :
        radius = radius_orig
        
    tbl_coord, tbl_flux = simbad_target_name_resolver(path, center)
    lambda_list = []
    dlambda_list = []
    
    
    if name == '0' :
        name = "III/284/allstars"
        url = f"https://vizier.cds.unistra.fr/viz-bin/votable?-source={name}&-c={center}&-c.rm={radius}&-out.all=1"
        system_list = ["2MASS", "2MASS", "2MASS", "Washington", "Washington", "Spitzer/IRAC", "Spitzer/IRAC", "Spitzer/IRAC", "Spitzer/IRAC", "WISE"]
        filter_list = ["J","H","Ks","M","T2","3.6","4.5","5.8","8.0","W2"]
        column_list = ["Jmag","Hmag","Ksmag","Mmag","T2mag","_3.6mag","_4.5mag","_5.8mag","_8.0mag","_4.5magW"]
        
        table = Table.read(url, format="votable")
        rac_list = table['RAJ2000']
        deg_list = table['DEJ2000']
        ID_list = table['ID']
        tbl_ = Table([ID_list, rac_list, deg_list], names=("ID","Ra","De"))
        if filtre_SED == "all" :
            nb_filter = len(filter_list)
        else :
            print("\nAvailable filters :\n0 : Jmag\n1 : Hmag\n2 : Ksmag\n3 : Mmag\n4 : T2mag\n5 : 3.6mag\n6 : 4.5mag\n7 : 5.8mag\n8 : 8.0mag\n9 : 4.5magW \n")
            nb_filter = int(input("[Number of filters you want (max 10)] "))        
        if nb_filter >= len(filter_list) :
            nb_filter = len(filter_list)
            for i in range(nb_filter) :
                tbl_ = all_filter(i, lambda_list, dlambda_list, filter_list, column_list, system_list, table, tbl_)
        else :
            filter_choose = []
            iteration = 1
            while iteration < nb_filter + 0.1 :
                tbl_, filter_choose = while_filter(lambda_list, dlambda_list, filter_list, column_list, system_list, table, tbl_, filter_choose, iteration)
                iteration += 1
        
    
    
    elif name == '1' :
        
        name = "II/340/xmmom2_1"
        url = f"https://vizier.cds.unistra.fr/viz-bin/votable?-source={name}&-c={center}&-c.rm={radius}&-out.all=1"
        system_list = ["XMM-OT"]
        filter_list = ["V","V","V","B","V"]
        column_list = ["UVM2mag","UVW1mag","Umag","Bmag","Vmag"]
        
        table = Table.read(url, format="votable")
        rac_list = table['RAICRS']
        deg_list = table['DEICRS']
        ID_list = table['ID']
        tbl_ = Table([ID_list, rac_list, deg_list], names=("ID","Ra","De"))
            
        if filtre_SED == "all" :
            nb_filter = len(filter_list)
        else :
            print("\nAvailable filters :\n0 : UVM2mag\n1 : UVW1mag\n2 : Umag\n3 : Bmag\n4 : Vmag\n")
            nb_filter = int(input("[Number of filters you want (max 5)] "))
        
        if nb_filter >= len(filter_list) :
            nb_filter = len(filter_list)
            for i in range(nb_filter) : 
                tbl_ = all_filter(i, lambda_list, dlambda_list, filter_list, column_list, system_list, table, tbl_)
        else :
            filter_choose = []
            iteration = 1
            while iteration < nb_filter + 0.1 :
                tbl_, filter_choose = while_filter(lambda_list, dlambda_list, filter_list, column_list, system_list, table, tbl_, filter_choose, iteration)
                iteration += 1


    elif name == '2' :
        name = "II/262/batc"
        url = f"https://vizier.cds.unistra.fr/viz-bin/votable?-source={name}&-c={center}&-c.rm={radius}&-out.all=1"
        system_list = ["BATC"]
        filter_list = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "m", "n", "o", "p"]
        column_list = ["aMag", "bMag", "cMag", "dMag", "eMag", "fMag", "gMag", "hMag", "iMag", "jMag", "kMag", "mMag", "nMag", "oMag", "pMag"]
        
        table = Table.read(url, format="votable")
        
        rac_list, deg_list = [], []
        for i in range(len(table['RAJ2000'])) :
            raj = ((float(table['RAJ2000'][i][0:2]) + float(table['RAJ2000'][i][3:5])/60 + float(table['RAJ2000'][i][6:])/3600)/24) * 360
            dej = float(table['RAJ2000'][i][0:3]) + float(table['RAJ2000'][i][4:6])/60 + float(table['RAJ2000'][i][7:])/3600
            rac_list.append(round(raj,8))
            deg_list.append(round(dej,8))
        
        ID_list = table['Seq']
        tbl_ = Table([ID_list, rac_list, deg_list], names=("ID","Ra","De"))
        
        if filtre_SED == "all" :
            nb_filter = len(filter_list)
        else :
            print("\nAvailable filters :\n0 : aMag\n1 : bMag\n2 : cMag\n3 : dMag\n4 : eMag\n5 : fMag\n6 : gMag\n7 : hMag\n8 : iMag\n9 : jMag\n10 : kMag\n11 : mMag\n12 : nMag\n13 : oMag\n14 : pMag\n")
            nb_filter = int(input("[Number of filters you want (max 15)] "))        
        if nb_filter >= len(filter_list) :
            nb_filter = len(filter_list)
            for i in range(nb_filter) :
                tbl_ = all_filter(i, lambda_list, dlambda_list, filter_list, column_list, system_list, table, tbl_)
        else :
            filter_choose = []
            iteration = 1
            while iteration < nb_filter + 0.1 :
                tbl_, filter_choose = while_filter(lambda_list, dlambda_list, filter_list, column_list, system_list, table, tbl_, filter_choose, iteration)
                iteration += 1
                

    else :
        print("[The name is not register, you will have to give all the collumn name]")
        
        url = f"https://vizier.cds.unistra.fr/viz-bin/votable?-source={name}&-c={center}&-c.rm={radius}&-out.all=1"
        
        try :
            table = Table.read(url, format="votable")
        except :
            print("\n[No table found, one of the name is wrong.]")
            tbl_ = Table([[nan], [nan], [nan]], names=("ID","Ra","De"))
            nb_filter = 0
            lambda_list = []
            radius = 0
            return tbl_, tbl_coord, nb_filter, lambda_list, radius

        name_ra = input("[Name of the Radial Ascention column] ")
        name_de = input("[Name of the Declination column] ")
        name_id = input("[Name of the stars' ID column] ")
        
        try :
            rac_list = table[name_ra]
            deg_list = table[name_de]
            ID_list = table[name_id]
        except :
            print("\n[Impossible to read the data, The name Radial Ascention or Declination or stars' ID is wrong.]")
            tbl_ = Table([[nan], [nan], [nan]], names=("ID","Ra","De"))
            nb_filter = 0
            lambda_list = []
            radius = 0
            return tbl_, tbl_coord, nb_filter, lambda_list, radius
        
        tbl_ = Table([ID_list, rac_list, deg_list], 
                          names=("ID","Ra","De"))
        
        Continue = "y"
        while Continue == "y" :            
            lambda_list.append([])
            dlambda_list.append([])
            column = int(input(f"[Column name in Vizier] "))
            try :
                column_err = "e_"+column
                try_table = table[column_err]
            except :
                column_err = "e"+column
                
            try :
                table[column]
                table[column_err]
            except :
                print(f"[ERROR : Choose a valid column name]")
                column = int(input(f"[Column name in Vizier] "))
                
            filter_ = input("[filter corresponding to this column]")
            system = input("[system use for this filter]")

            lambda_, dlambda, Fmag = search_vega_filter(system, filter_)
            flux_jsky_list = []
            err_flux_jsky_list = []

        
            for i in range(len(table[column])) :
                table[column].fill_value = nan
                table[column_err].fill_value = nan
                     
                mag_err = table[column_err].filled()[i]
                mag = table[column].filled()[i]
                        
                mag_neg = mag - mag_err
                mag_pos = mag + mag_err
                flux = to_jsky(Fmag,mag)
                error_neg = round(flux - to_jsky(Fmag, mag_neg), 4)
                error_pos = round(flux - to_jsky(Fmag, mag_pos), 4)
                flux_error = f"+{error_pos}, {error_neg}"
                        
                flux_jsky_list.append(flux)
                err_flux_jsky_list.append(flux_error)
                       
                lambda_list[iteration-1].append(lambda_)
                dlambda_list[iteration-1].append(dlambda)
                        
            tbl_.add_column(flux_jsky_list*unit_jsky, name=f"{column}")
            tbl_.add_column(err_flux_jsky_list, name=f"{column_err} (jsky)")
                    
            for j in range(len(flux_jsky_list)) :
                if np.isnan(flux_jsky_list[j]) == True :
                    print(f"[WARNING : There are some empty data in the filter {filter_}.]\n")
                    break
    
    tbl_.write(f'{path}tbl_flux_vizier.ecsv', overwrite=True)           
    
    return tbl_, tbl_coord, nb_filter, lambda_list, dlambda_list, radius, center


########################################################################################################################################################################################################
            # SIMBAD AND VIZIER PHOTOMETRIC RESOLVER :
########################################################################################################################################################################################################


def plot_errorbar (nb_filter, tbl_, lambda_list, dlambda_list) :
    for j in range(nb_filter) :
        if np.nansum(tbl_[tbl_.keys()[2*j+3]]) > 0 :
            print(f"[Lambda for the filter {tbl_.keys()[2*j+3]} is {lambda_list[j][0]}]")
            yerr_up, yerr_down = [], []
            yerr_array = np.array(tbl_[tbl_.keys()[2*j+4]])
            for k in range(len(yerr_array)):    
                yerr_up.append(float( yerr_array[k].split(",")[0] ))
                yerr_down.append(float( yerr_array[k].split(",")[1] ))                 
                
            plt.errorbar(lambda_list[j], tbl_[tbl_.keys()[2*j+3]], xerr=dlambda_list[j], yerr=(yerr_down,yerr_up), barsabove=True, capsize=3, label=f' {tbl_.keys()[2*j+3]}')
        else :
            print(f"[WARNING : There is no value for {tbl_.keys()[2*j+3]}.]")

def plot(nb_filter, tbl_, lambda_list) :
    for j in range(nb_filter) :
        if np.nansum(tbl_[tbl_.keys()[2*j+3]]) > 0 :
            print(f"[Lambda for the filter {tbl_.keys()[2*j+3]} is {lambda_list[j][0]}]")
            plt.scatter(lambda_list[j], tbl_[tbl_.keys()[2*j+3]])
        else :
            print(f"[WARNING : There is no value for {tbl_.keys()[2*j+3]}.]")

def SED(path, nb_catalogue=0) :
    list_cata = ["III/284/allstars","II/340/xmmom2_1","II/262/batc"]
    error_bar = input("[Do you want to display the error bar in the SED plot ? (y/n)] ")
    while error_bar != "y" and error_bar != "n" :
        print("[Please enter a valid answer]")
        error_bar = input("[Do you want to display the error bar in the SED plot ? (y/n)] ")

    center_name = input("[Conesearch's center name (e.g. M77, HD1, NGC2264)] ")
    radius_conesearch = input("[Conesearch's radius (arcmin)] ")
    
    print("\n\nCONESEARCH : \n\n")
    plt.figure()
    
    for i in range(0,nb_catalogue):
        print("\nCatalogue {list_cata[i]}\n")
        tbl_, coord_origin, nb_filter, lambda_list, dlambda_list, radius_orig, center_name = vizier_cone_search(path, "all", i, center_name, radius_conesearch)
        plt.scatter(tbl_["Ra"], tbl_["De"], s=5)
        
    rac_origin, deg_origin = coord_origin["Ra "], coord_origin["Dec "]
    plt.scatter(rac_origin, deg_origin, c='red', label=f"{center_name}")

    plt.xlabel("Right Ascension (degres)")
    plt.ylabel("Declinaison (degres)")
    radius_conesearch = float(radius_conesearch)
    plt.title(f"Conesearch {radius_conesearch}arcmin \n around {center_name}", fontsize=20)
    plt.legend()
    
    plt.savefig(f"{path}Conesearch.png")
    
    print("\n\nSED : \n\n")
    plt.figure()
    
    if error_bar == "y" :                
        for i in range(0,nb_catalogue):
            print(f"\n[Catalogue {list_cata[i]}]\n")
            tbl_, coord_origin, nb_filter, lambda_list, dlambda_list, radius, center_name = vizier_cone_search(path, "all", i, center_name, radius_orig)
            plot_errorbar(nb_filter, tbl_, lambda_list, dlambda_list)
    else :
        for i in range(0,nb_catalogue):
            print(f"\n[Catalogue {list_cata[i]}]\n")
            tbl_, coord_origin, nb_filter, lambda_list, dlambda_list, radius, center_name = vizier_cone_search(path, "all", i, center_name, radius_orig)
            plot(nb_filter, tbl_, lambda_list)
    
    plt.yscale("log")
    plt.xlabel("lambda (micrometer)")
    plt.ylabel("Flux (Jsky)")
    plt.title("Spectral Distribution of Energy", fontsize=20)
    
    plt.savefig(f"{path}SED.png")
    
    plt.show()

