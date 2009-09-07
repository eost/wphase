/* Convert date to julian days */
int 
yyyymmdd2jjj(int year, int month, int day)
{

  int ndays[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
  int k, jul;
      
  if ( year%4   == 0 ) 
	ndays[1] += 1;
  if ( year%100 == 0 ) 
	ndays[1] -= 1;
  if ( year%400 == 0 ) 
	ndays[1] += 1;
      
  jul = day;
  for(k=0; k<month-1; k++)
	jul += ndays[k];
  
  return(jul);
}
