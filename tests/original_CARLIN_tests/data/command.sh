for x in *.gz
 do
   echo $x
   wc -l $x
done
