#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 09:57:06 2022

@author: rvancoellie
"""

####################################################################################################
            # IMPORT :
####################################################################################################

import urllib.request as http
from astropy.table import Table, Column
import numpy as np
import matplotlib.pyplot as plt
from cython_code import c_to_jsky as to_jsky
from cython_code import c_search_vega_filter as search_vega_filter

########################################################################################################################################################################################################
            # SIMBAD TARGET RESOLVER :
########################################################################################################################################################################################################

def Diplay_simbad_target_name_resolver(nameornot, name="Giacomo") :
    tbl_coordinate, tbl_flux = simbad_target_name_resolver(nameornot, name)
    Ra_deg = tbl_coordinate[0]
    Dec_deg = tbl_coordinate[1]
    a=input("avant coord")
    print("\nCoordinate : ")
    print("Ra : ", Ra_deg, "degres, \nDec : ", Dec_deg, "degres.")
    a=input("après coord")

    a=input("avant tbl flux")
    print(tbl_flux)
    a=input("après tbl flux")
    return 1

def simbad_target_name_resolver(nameornot, name="Giacomo") :
    if nameornot == 0 :
        TARGET = input('Target name : ')
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
    

    with http.urlopen(url) as fd :
        for line in fd:
            if "Coordinates(ICRS,ep=J2000,eq=2000)" in line.decode('UTF-8'):
                r_line = line.decode('utf-8').strip().split()
                Ra_deg = ((float(r_line[1])  + float(r_line[2])/60 + float(r_line[3])/3600)/24) * 360
                Dec_deg = float(r_line[4])  + (float(r_line[5])/60 + float(r_line[6])/3600)
                tbl_coordinate = [Ra_deg, Dec_deg]
                
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
                        
                        error_neg = round(flux - flux_neg, 8)
                        error_pos = round(flux - flux_pos, 8)
                        
                        flux_error = f"+{error_pos}, {error_neg}"
                        
                    else :
                        print(f"There is no flux error in Simbad data for the filter {filter_name}")
                        flux_error = f"+{nan}, -{nan}"
                    
                    ar_filter_name.append(filter_name)
                    ar_flux.append(flux)
                    ar_flux_error.append(flux_error)
                    ar_lambda.append(lambda_)
                    ar_lambda_error.append(dlambda)

                    
    tbl_flux = Table([ar_filter_name, ar_flux, ar_flux_error, ar_lambda ,ar_lambda_error], 
                  names=("filter","flux (jsky)","flux error (jsky)","Lambda","Lambda error"))
    
    return tbl_coordinate, tbl_flux

########################################################################################################################################################################################################
            # VIZIER CONE SEARCH :
########################################################################################################################################################################################################

def vizier_cone_search() :
    nan = float("NaN")
    name = input('name : ')
    center = input('center : ')
    radius = input('radius (arcmin) : ')
    url = f"https://vizier.cds.unistra.fr/viz-bin/votable?-source={name}&-c={center}&-c.rm={radius}&-out.all=1"
    
    tbl_coord, tbl_flux = simbad_target_name_resolver(1,center)
    lambda_list = []
    
    if name == 'III/284/allstars' :
        
        table = Table.read(url, format="votable")
        rac_list = table['RAJ2000']
        deg_list = table['DEJ2000']
        ID_list = table['ID']
        tbl_ = Table([ID_list, rac_list, deg_list], 
                          names=("ID","Ra","De"))
        
        nb_filter = int(input("How many filter do you want to select : "))
        i=0
        for i in range(nb_filter) :
            lambda_list.append([])
            
        iteration = 1
        print("\n0 : Jmag\n1 : Hmag\n2 : Ksmag\n3 : Mmag\n4 : T2mag\n5 : 3.6mag\n6 : 4.5mag\n7 : 5.8mag\n8 : 8.0mag\n9 : 4.5magW \n")
        
        while iteration < nb_filter + 0.1 :
            System_list = ["2MASS", "2MASS", "2MASS", "Washington", "Washington", "Spitzer/IRAC", "Spitzer/IRAC", "Spitzer/IRAC", "Spitzer/IRAC", "WISE"]
            filter_list = ["J","H","Ks","M","T2","3.6","4.5","5.8","8.0","W2"]
            column_list = ["Jmag","Hmag","Ksmag","Mmag","T2mag","_3.6mag","_4.5mag","_5.8mag","_8.0mag","_4.5magW"]
            column_nb = int(input(f"Choose the filter number {round(iteration,1)} you want : "))
            system = System_list[column_nb]
            filter_ = filter_list[column_nb]
            column = column_list[column_nb]
            print(f"for the column {column} the system use is {system} with the filter named {filter_}.")

            lambda_, dlambda, Fmag = search_vega_filter(system, filter_)
            flux_jsky_list = []


            for i in range(len(table[column])) :
                table[column].fill_value = nan
                mag = table[column].filled()[i]
                flux_jsky_list.append(to_jsky(Fmag,mag))
                lambda_list[iteration-1].append(lambda_)


            tbl_.add_column(flux_jsky_list, name=f"Flux {filter_} (Jsky)")
            
            for j in range(len(flux_jsky_list)) :
                if np.isnan(flux_jsky_list[j]) == True :
                    print(f"There are some empty data in the filter {filter_} you choose.\n")
                    break
            
            
            iteration += 1
    
    elif name == 'II/340/xmmom2_1' :
        
                
        table = Table.read(url, format="votable")
        rac_list = table['RAICRS']
        deg_list = table['DEICRS']
        ID_list = table['ID']        
        tbl_ = Table([ID_list, rac_list, deg_list], 
                          names=("ID","Ra","De"))
        
        nb_filter = int(input("How many filter do you want to select (max 5): "))
        i=0
        for i in range(nb_filter) :
            lambda_list.append([])
            
        iteration = 1
        print("\n0 : UVM2mag\n1 : UVW1mag\n2 : Umag\n3 : Bmag\n4 : Vmag\n")
        
        while iteration < nb_filter + 0.1 :
            System = "XMM-OT"
            filter_list = ["V","V","V","B","V"]
            column_list = ["UVM2mag","UVW1mag","Umag","Bmag","Vmag"]
            column_nb = int(input(f"Choose the filter number {round(iteration,1)} you want : "))
            filter_ = filter_list[column_nb]
            column = column_list[column_nb]
            print(f"for the column {column} the system use is {System} with the filter named {filter_}.")

            lambda_, dlambda, Fmag = search_vega_filter(System, filter_)
            table = Table.read(url, format="votable")
            flux_jsky_list = []

            for i in range(len(table[column])) :
                table[column].fill_value = nan
                mag = table[column].filled()[i]
                flux_jsky_list.append(to_jsky(Fmag,mag))
                lambda_list[iteration-1].append(lambda_)


            tbl_.add_column(flux_jsky_list, name=f"{column} (Jsky)")

            for j in range(len(flux_jsky_list)) :
                if np.isnan(flux_jsky_list[j]) == True :
                    print(f"There are some empty data in the filter {filter_} you choose.\n")
                    break
            
            iteration += 1

    else :
        print("Name not known, we can't resolve the system")
        tbl_ = nan
                
    return tbl_, tbl_coord, nb_filter, lambda_list

####################################################################################################
            # SIMBAD AND VIZIER PHOTOMETRIC RESOLVER :
####################################################################################################

def photo() :
    tbl_, coord_origin, nb_filter, lambda_list = vizier_cone_search()
    
    rac_origin, deg_origin = coord_origin[0], coord_origin[1]
    
    i=0
    for i in range(len(tbl_["Ra"])) :
        if tbl_["Ra"][i] > 180 :
            tbl_["Ra"][i] -= 360
    
    if rac_origin > 180 :
        rac_origin -= 360
    
    plt.figure()
    
    plt.scatter(rac_origin, deg_origin, c='red')
    plt.scatter(tbl_["Ra"], tbl_["De"], c='blue')

    plt.xlabel("Right Ascension (degres)")
    plt.ylabel("Declinaison (degres)")
    
    plt.show()
    
    
    plt.figure()
    
    i=0
    for i in range(nb_filter) :
        print(f"Lambda for the filter {tbl_.keys()[i+3]} is {lambda_list[i][0]}")
        plt.scatter(lambda_list[i], tbl_[tbl_.keys()[i+3]])

    plt.xlabel("lambda (?)")
    plt.ylabel("Flux (Jsky)")

    plt.show()

