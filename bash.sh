st_year=$1
ed_year=$2
more=$3
echo "end year = $st_year"
echo "star year = $ed_year"
echo "more = $more"
path1="/raid61/WRF_reinit_GLDASsnow_SKINsst/Output"
path2="/raid61/WRF_reinit_GLDASsnow_SKINsst/Output"
echo "path1 = $path1"
echo "path2 = $path2"
python prcp.py $st_year $ed_year $more $path1 $path2
python T2.py $st_year $ed_year $more $path1 $path2
python u10.py $st_year $ed_year $more $path1 $path2
python v10.py $st_year $ed_year $more $path1 $path2
python Q2.py $st_year $ed_year $more $path1 $path2

var_names=("PBLH" "ALBEDO" "HFX" "QFX" "PSFC" "LH" "GLW" "SWDOWN" "SNOWH" "SNOW" "GRDFLX")
units=("m" "-" "W*m-2" "kg*m-2*s-1" "Pa" "W*m-2" "W*m-2" "W*m-2" "m" "kg*m-2" "W*m-2")
file_var_names=("pblh" "albedo" "hfx" "qfx" "psfc" "lh" "glw" "swdown" "snowh" "snow" "grdflx")
for ((i=0;i<${#var_names[@]};i++))
do
    echo ${var_names[$i]} ${units[$i]} ${file_var_names[$i]}
    python surface_instant_var.py $st_year $ed_year $more $path1 $path2 ${var_names[$i]} ${units[$i]} ${file_var_names[$i]}
done


