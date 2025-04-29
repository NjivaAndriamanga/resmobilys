BEGIN { FS=OFS="\t" }
NR==1 {
    print ; 
}

/Candidate/ {print}
/^Predicted/ {flag=1} /Candidate/ {flag=0} flag && !/^Candidate/ {print}
