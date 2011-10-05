#
#                            Fuzzy String Match
#
#   Copyright 2010 Kiyoka Nishiyama
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
module FuzzyStringMatch

  class JaroWinkler
    def create( type = :pure )     # factory method
      case type
      when :native
        JaroWinklerNative.new
      when :pure
        JaroWinklerPure.new
      end
    end
  end
  
  class Bitap
    def create(type = :native)
      case type
      when :native
        BitapNative.new
      when :pure
        BitapPure.new
      end
    end
  end

  class Levenshtein
    def create(type = :native)
      case type
      when :native
          LevenshteinNative.new
      when :pure
        LevenshteinPure.new
      end

    end
  end
  
  class Ngram
    def create( type = :native )
      case type
      when :native
        NgramNative.new
      when :pure
        NgramPure.new
      end
    end
  end
  
  class JaroWinklerPure
    THRESHOLD = 0.7

    def getDistance( s1, s2 )
      a1 = s1.split( // )
      a2 = s2.split( // )

      if s1.size > s2.size
        (max,min) = a1,a2
      else
        (max,min) = a2,a1
      end

      range = [ (max.size / 2 - 1), 0 ].max
      indexes = Array.new( min.size, -1 )
      flags   = Array.new( max.size, false )

      matches = 0;
      (0 ... min.size).each { |mi|
        c1 = min[mi]
        xi = [mi - range, 0].max
        xn = [mi + range + 1, max.size].min

        (xi ... xn).each { |i|
          if (not flags[i]) && ( c1 == max[i] )
            indexes[mi] = i
            flags[i] = true
            matches += 1
            break
          end
        }
      }

      ms1 = Array.new( matches, nil )
      ms2 = Array.new( matches, nil )

      si = 0
      (0 ... min.size).each { |i|
        if (indexes[i] != -1)
          ms1[si] = min[i]
          si += 1
        end
      }

      si = 0
      (0 ... max.size).each { |i|
        if flags[i]
          ms2[si] = max[i]
          si += 1
        end
      }

      transpositions = 0
      (0 ... ms1.size).each { |mi|
        if ms1[mi] != ms2[mi]
          transpositions += 1
        end
      }

      prefix = 0
      (0 ... min.size).each { |mi|
        if s1[mi] == s2[mi]
          prefix += 1
        else
          break
        end
      }

      if 0 == matches
        0.0
      else
        m = matches.to_f
        t = (transpositions/ 2)
        j = ((m / s1.size) + (m / s2.size) + ((m - t) / m)) / 3.0;
        return j < THRESHOLD ? j : j + [0.1, 1.0 / max.size].min * prefix * (1 - j)
      end
    end
  end

  class BitapPure
  end

  class LevenshteinPure
  end
  
  require 'inline'
  class JaroWinklerNative
    inline do |builder|
      builder.include '<iostream>'
      builder.add_compile_flags '-x c++', '-lstdc++'
      builder.c_raw 'int max( int a, int b ) { return ((a)>(b)?(a):(b)); }'
      builder.c_raw 'int min( int a, int b ) { return ((a)<(b)?(a):(b)); }'
      builder.c_raw 'double double_min( double a, double b ) { return ((a)<(b)?(a):(b)); }'
      builder.c '
double getDistance( char *s1, char *s2 )
{
  char *_max;
  char *_min;
  int _max_length = 0;
  int _min_length = 0;
  if ( strlen(s1) > strlen(s2) ) {
    _max = s1; _max_length = strlen(s1);
    _min = s2; _min_length = strlen(s2);
  }
  else {
    _max = s2; _max_length = strlen(s2);
    _min = s1; _min_length = strlen(s1);
  }
  int range = max( _max_length / 2 - 1, 0 );

  int indexes[_min_length];
  for( int i = 0 ; i < _min_length ; i++ ) {
    indexes[i] = -1;
  }

  int flags[_max_length];
  for( int i = 0 ; i < _max_length ; i++ ) {
    flags[i] = 0;
  }
  int matches = 0;
  for (int mi = 0; mi < _min_length; mi++) {
    char c1 = _min[mi];
    for (int xi = max(mi - range, 0), xn = min(mi + range + 1, _max_length); xi < xn; xi++ ) {
      if (!flags[xi] && (c1 == _max[xi])) {
	indexes[mi] = xi;
	flags[xi] = 1;
	matches++;
	break;
      }
    }
  }

  char ms1[matches];
  char ms2[matches];
  int ms1_length = matches;

  for (int i = 0, si = 0; i < _min_length; i++) {
    if (indexes[i] != -1) {
      ms1[si] = _min[i];
      si++;
    }
  }
  for (int i = 0, si = 0; i < _max_length; i++) {
    if (flags[i]) {
      ms2[si] = _max[i];
      si++;
    }
  }
  int transpositions = 0;
  for (int mi = 0; mi < ms1_length; mi++) {
    if (ms1[mi] != ms2[mi]) {
      transpositions++;
    }
  }
  int prefix = 0;
  for (int mi = 0; mi < _min_length; mi++) {
    if (s1[mi] == s2[mi]) {
      prefix++;
    } else {
      break;
    }
  }

  double m = (double) matches;
  if (matches == 0) {
    return 0.0;
  }
  int t = transpositions / 2;
  double j = ((m / strlen(s1) + m / strlen(s2) + (m - t) / m)) / 3;
  double jw = j < 0.7 ? j : j + double_min(0.1, 1.0 / _max_length) * prefix
    * (1 - j);
  return jw;
}'
    end
  end

  class BitapNative
    inline :C do |builder|
      builder.include '<stdlib.h>'
      builder.include '<string.h>'
      builder.include '<limits.h>'
      builder.c "const char *search(const char *text, const char *pattern, int k)
       {
           const char *failure = 'Failed';
           const char *result = NULL;
           int m = strlen(pattern);
           unsigned long *R;
           unsigned long pattern_mask[CHAR_MAX+1];
           int i, d;
           if (pattern[0] == '\0') return text;
           if (m > 31) return \"The pattern is too long!\";
           /* Initialize the bit array R */
           R = malloc((k+1) * sizeof *R);
           for (i=0; i <= k; ++i)
               R[i] = ~1;
           /* Initialize the pattern bitmasks */
           for (i=0; i <= CHAR_MAX; ++i)
               pattern_mask[i] = ~0;
           for (i=0; i < m; ++i)
               pattern_mask[pattern[i]] &= ~(1UL << i);

           for (i=0; text[i] != '\0'; ++i) {
               /* Update the bit arrays */
               unsigned long old_Rd1 = R[0];

               R[0] |= pattern_mask[text[i]];
               R[0] <<= 1;

               for (d=1; d <= k; ++d) {
                   unsigned long tmp = R[d];
                   /* Substitution is all we care about */
                   R[d] = (old_Rd1 & (R[d] | pattern_mask[text[i]])) << 1;
                   old_Rd1 = tmp;
               }

               if (0 == (R[k] & (1UL << m))) {
                   result = (text+i - m) + 1;
                   break;
               }
           }
           free(R);
           if(result == NULL)
             result = \"0\";
           else
             result = \"1\";
           return result;
       }"
    end
  end
  
  class LevenshteinNative
    inline :C do |build|
      build.include '<stdlib.h>'
      build.include '<string.h>'
      build.c '
      int search(char *s,char*t);
      int minimum(int a,int b,int c);

      int search(char *s,char*t)
      /*Compute levenshtein distance between s and t*/
      {
        //Step 1
        int k,i,j,n,m,cost,*d,distance;
        n=strlen(s); 
        m=strlen(t);
        if(n!=0&&m!=0)
        {
          d=malloc((sizeof(int))*(m+1)*(n+1));
          m++;
          n++;
          //Step 2	
          for(k=0;k<n;k++)
      	d[k]=k;
          for(k=0;k<m;k++)
            d[k*n]=k;
          //Step 3 and 4	
          for(i=1;i<n;i++)
            for(j=1;j<m;j++)
      	{
              //Step 5
              if(s[i-1]==t[j-1])
                cost=0;
              else
                cost=1;
              //Step 6			 
              d[j*n+i]=minimum(d[(j-1)*n+i]+1,d[j*n+i-1]+1,d[(j-1)*n+i-1]+cost);
            }
          distance=d[n*m-1];
          free(d);
          return distance;
        }
        else 
          return -1; //a negative return value means that one or both strings are empty.
      }

      int minimum(int a,int b,int c)
      /*Gets the minimum of three values*/
      {
        int min=a;
        if(b<min)
          min=b;
        if(c<min)
          min=c;
        return min;
      }'
    end
  end
  
  class NgramNative
    inline :C do |builder|
      #builder.include '<ctypes.h>'
      builder.include '<stdlib.h>'
      builder.include '<string.h>'
      builder.include '<math.h>'
      builder.c_raw '
      #include <stdio.h>
      #include <math.h>
      #include <string.h>
      #include <strings.h>
      #include <stdlib.h>

      typedef struct{
      	char *theWord;
      	int theCount;
      } GRAM;

      int length1 = -1;
      int length2 = -1;
      GRAM **grams1 = NULL;
      GRAM **grams2 = NULL;

      //int common(result *one, result *two);
      //int matches(result *one, result *two);
      //int contains(const char *c, const result *results);
      //int processString(const char *c, const int n, result *res);
      //int my_substr(int from, int to, char* str, char* substr);
      //double getSimilarity(const char *termOne, const char *termTwo, int n);

      int addGram1(GRAM g){
      	if(length1 == -1){
      		length1 = 0;
      	}else{
      		length1 += 1;
      	}

      	void *_tmp = realloc(grams1, (length1 * sizeof(GRAM)));
      	if (!_tmp)
      	{
      		fprintf(stderr, "ERROR: Couldn't realloc memory!\n");
      		return(-1);
      	}	grams1 = (GRAM*)_tmp;
      	*grams1[length1] = g;
      	return length1;
      }//addGram1

      int addGram2(GRAM g){
      	if(length2 == -1){
      		length2 = 0;
      	}else{
      		length2 += 1;
      	}

      	void *_tmp = realloc(grams2, (length2 * sizeof(GRAM)));
      	if (!_tmp)
      	{
      		fprintf(stderr, "ERROR: Couldn't realloc memory!\n");
      		return(-1);
      	}	grams2 = (GRAM*)_tmp;
      	int s1 = length2;
      	*grams2[length2] = g;
      	char *s = *grams2[length2]->theWord;
      	return length2;
      }//addGram2

      int common(GRAM one[], GRAM two[]){
      	//char *s = one[1].theWord;
      	int res = 0, i, j;
      	char *word = grams2[1]->theWord;
      	for(i = 0; i <= length1; i++){
      		for(j = 0; j <= length2; j++){
      			char *s = two[j].theWord;
      			char *m = one[i].theWord;
      			if(!strcasecmp(one[i].theWord, two[j].theWord)){
      				res++;
      			}
      		}
      	}
      	return res;
      	//return 0;
      }//end common

      int matches(GRAM one[], GRAM two[]){
      	int i, j;
      	for(i = 0; i <= length2; i++){
      		int index = -1;
      		int found = 0;
      		for(j = 0; j < length1 && !found; j++){
      			if(!strcasecmp(two[i].theWord, one[j].theWord)){
      				found = 1;
      			}
      			index = j;
      		}

      		if(!found){
      			one[length1 - 1] = two[index];
      		}
      	}

      	return length1;
      }//end matches

      int contains(const char *c, const GRAM results[]){
      	int i;
      	int size = sizeof(results)/sizeof(results[0]);
      	for(i = 0; i < size; i++){
      		GRAM s = results[i];
      		if(s.theWord != '\0'){
      			if(!strcasecmp(results[i].theWord, c)){
      				return 1;
      			 }
      		}
      	}
      	return 0;
      }//end contains

      int processString(const char *c, const int n, int flag){
      //	int count = ceil((sizeof(c) / n));
      	int index = 0;
      	int i, j, k;
      //	result r[count];
      //	res = &r;

      	for(i = 0; i < strlen(c) - 1; i++){
      		char tmp[n];
      //		char final[n];
      		for(j = i, k = 0; j < (i+n); j++){
      			tmp[k] = c[j];
      			k++;
      		}
      //		final = tmp;
      		if(flag == 1 && !contains(tmp, grams1)){
      			//GRAM st_tmp = {"Jo", length1};
      			GRAM st_tmp = {tmp, length1};
      			length1 = addGram1(st_tmp);	
      			//res[index] = &st_tmp;
      			//res[index]->theWord = "Jo";
      			//res[index]->theCount = 50;
      //			//int strlength = strlen(tmp);
      //			strncpy(*st_tmp.theWord, tmp, strlen(tmp));
      //			
      //			st_tmp.theCount = 23;
      //			result *res_tmp;
      //			res_tmp = &res[index];
      //			*res_tmp = st_tmp;
      		}
      		if(flag == 2 && !contains(tmp, grams2)){
      			//GRAM st_tmp = {"Jo", length2};
      			GRAM st_tmp = {tmp, length2};
      			length2 = addGram2(st_tmp);	
      		}
      	}//for
      	return 1;
      }//end processString

      int my_substr(int from, int to, char* str, char* substr){
      	int i, k;

      	if(to > from){
      		return 1;
      	}
      	int substr_size = (to - from ) + 2;
      	substr[substr_size];
      	for(i = from, k = 0; i <= to && k < substr_size; i++){
      		substr[k] = str[i];
      		k++;
      	}  
      	return 1;
      }//end my_substr

      const double *getSimilarity(const char *termOne, const char *termTwo, int n){
      	int count_1 = ceil((strlen(termOne) - n));
      	int count_2 = ceil((strlen(termTwo) - n));
      //	int myInt[4];
      //	//result *res1[count_1];
      //	result **res1;
      //	res1 = malloc(((sizeof(char) * n) + sizeof(int)) * count_1);
      //       //result* res1 = malloc(count_1 * (sizeof(int) + sizeof(result) * n));
      //	//res1 = malloc(count_1 * sizeof(result));
      //	result **res2;
      //	res2 = malloc(((sizeof(char) * n) + sizeof(int)) * count_2);
      	processString(termOne, n, 1);
      	processString(termTwo, n, 2);
      	int c = common(grams1, grams2);
      	int m = matches(grams1, grams2);
      	double *sim;
      	*sim = (double) c / (double) m;
      	const double *val = sim;
      	return val;
      }

      int main (int argc, const char * argv[]) {
          // insert code here...
          printf("Hello, World!\n");
      	double *s = getSimilarity("John", "Johnny", 2);
      	printf("The Double Is %f", *s);
          return 0;
      }

      '
    end
  end
end
