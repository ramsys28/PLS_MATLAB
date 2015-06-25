# PLS_MATLAB

Assume that for a single subject we have a NxN correlation matrix, called CM.

We should keep only the upper triangular part of the correlation matrix, using the following MATLAB command
CM_U = triu(CM,1);
Then we should reshape it in a single row using this command:
CM_U_res=reshape(CM_U,1,[]);
If we do this for all of our subjects then we will have a matrix with so many rows as our subjects. 
Lets call this matrix conns.
So conns contains in the n^{th} row the connexels of the n^{th}subjects.

Now we have to form the behavioral matrix. Each line of the behavioral matrix should contain the scores for a single subject.
This means that in the first line of the matrix we will have the behavioral scores of the first subject. 
In the second line the same for the second subject, so on...
Let's call this matrix behav

ATTENTION: THE CONNS AND BEHAV MATRICES SHOULD BE ALLIGNED ACCROSS THE SUBJECTS. THIS MEANS THAT THE FIRST LINE
SHOULD CONTAIN THE DATA OF THE FIRST SUBJECT (IN BOTH CONNS AND BEHAV MATRICES), THE SECOND LINE SHOULD CONTAIN THE DATA 
OF THE SECOND SUBJECTS AND SO ON.

After finishing with our matrices we just call the function
result=for_harvard(conns,behav)

and we have the results. 

