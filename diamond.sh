diamond makedb --in $1.faa -d $1.dmnd --quiet
diamond makedb --in $2.faa -d $2.dmnd --quiet

diamond blastp --query $1.faa -d $1.dmnd  -o $1.$1.m8 -f tab --min-score 50 --quiet
diamond blastp --query $1.faa -d $2.dmnd  -o $1.$2.m8 -f tab --min-score 50 --quiet
diamond blastp --query $2.faa -d $1.dmnd  -o $2.$1.m8 -f tab --min-score 50 --quiet
diamond blastp --query $2.faa -d $2.dmnd  -o $2.$2.m8 -f tab --min-score 50 --quiet
