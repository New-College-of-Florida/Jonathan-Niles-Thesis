RESOLUTIONS=("100kb" "200kb" "400kb" "800kb" "1000kb")
for i in "${RESOLUTIONS[@]}"; do ./script.sh $i; done;
