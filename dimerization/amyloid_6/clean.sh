#!/bin/bash

rm 1_water_min/mdinfo 1_water_min/min_wat.out 1_water_min/min_wat.rst
rm 2_water_eq/mdinfo  2_water_eq/mdcrd  2_water_eq/md_wat.out  2_water_eq/md_wat.rst  2_water_eq/min_wat.rst 
rm 3_sys_min/mdinfo 3_sys_min/md_wat.rst 3_sys_min/sys_min.out 3_sys_min/sys_min.rst
rm 4_sys_eq/heat.out 4_sys_eq/heat.rst 4_sys_eq/mdcrd 4_sys_eq/mdinfo 4_sys_eq/sys_min.rst
rm 5_meta_COM/COLVAR 5_meta_COM/heat.rst 5_meta_COM/HILLS 5_meta_COM/mdcrd 5_meta_COM/mdinfo  5_meta_COM/prod.out 5_meta_COM/test*
