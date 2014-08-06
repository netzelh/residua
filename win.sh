STAR=$1
PERIOD=$2
FREQ=$3

echo "read a $STAR.dat 3 1 2" > komendy.exec
echo "setp a $PERIOD" >> komendy.exec
echo 'fd a 5' >> komendy.exec
echo "dft a 0 10 10 win $FREQ" >> komendy.exec
echo "writefs a temp.txt" >> komendy.exec
echo 'exit' >> komendy.exec

fdecomp komendy.exec
