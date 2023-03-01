#!/bin/bash
dessin.gnu (){
        echo $1 $2 $3
    j=$(echo $1 | awk -F'.' '{print $1}')
    sed "s/input/$j/g" script_f.plt > script_fmod.plt
    sed -i -e "s/vmin/$2/g" script_fmod.plt
    sed -i -e "s/vmax/$3/g" script_fmod.plt
    gnuplot script_fmod.plt
}

#vmin=-3.1231550063348977
#vmax=3.1231550063348977
for i in $(ls vx*.plt); do
    dessin.gnu $i $vmin $vmax
done

#vmin=-3.1231550063348977
#vmax=3.1231550063348977
for i in $(ls vx-vx_a*.plt); do
    dessin.gnu $i $vmin $vmax
done

#vmin=-6.2522554097300338
#vmax=6.2522554097300338
for i in $(ls vz*.plt); do
    dessin.gnu $i $vmin $vmax
done

#vmin=-6.2522554097300338
#vmax=6.2522554097300338
for i in $(ls vz-vz_a*.plt); do
    dessin.gnu $i $vmin $vmax
done


#vmin=-13.480172263961662
#vmax=13.480172263961662
for i in $(ls pr*.plt); do
    dessin.gnu $i $vmin $vmax
done
    

#vmin=-12.204071419565713
#vmax=12.204071419565713
for i in $(ls pr_a*.plt); do
    dessin.gnu $i $vmin $vmax
done

#vmin=-11.52
#vmax=11.52
for i in $(ls pr-pr_a*.plt); do
    dessin.gnu $i $vmin $vmax
done


