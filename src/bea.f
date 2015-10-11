      subroutine rbea(n,m,a,istart,b,ib,ifin)
c     Reorder rows using BEA, bond energy algorithm.
      dimension a(n,m), b(n,m), ib(n), ifin(n)
c------------------------------------------------------------------------------
c     a(n,m)         input matrix, rows of which are to be permuted
c     istart         1st row to be placed
c     b(n,m)         permuted rows to be stored in this array
c     ib(n)          integer list giving permutation carried out
c     ifin(n)        book-keeping vector: is row still active, or has it been 
c                    placed (resp. = 1, = 0).
c------------------------------------------------------------------------------
c     Bond energy algorithm -- see:
c
c     (1) W.T. McCormick, P.J. Schweitzer and T.W. White, 
c         "Problem decomposition and data reorganization by a clustering 
c         technique", Oper. Res., vol. 20, pp. 993-1009, Sept./Oct. 1972.
c     (2) P. Arabie and L.J. Hubert, 
c         "The bond energy algorithm revisited", IEEE Trans. Syst. Man 
c         Cybern., vol. 20, pp. 268-274, 1990.
c     (3) P. Arabie, S. Schleutermann, J. Daws and L. Hubert, 
c         "Marketing applications of sequencing and partitioning of 
c         nonsymmetric and/or two-mode matrices", in W. Gaul and M. Schader, 
c         Eds., Data Analysis, Decision Support, and Expert Knowledge 
c         Representation in Marketing, Springer Verlag, 1988, pp. 215-224.
c  
c     Implemented by F. Murtagh, Sept. 1991.
c------------------------------------------------------------------------------
c
c     Flags to indicate if row already chosen; 1 = not yet chosen/placed.
      do 200 i = 1, n
 200     ifin(i) = 1
c
c     Anticipate 1st placement. 'nplace' = # rows placed. 'nrem' = # remaining.
      nplace = 1
      nrem   = n-1
c
c     Place 1st row
      do 300 j = 1, m
         b(1,j)     = a(istart,j)
         ifin(istart) = 0
         ib(nplace) = istart
 300  continue
c
 400  continue
      sim   = -100000.0
      insrt = 0
c     'nplace' rows have been placed.
c     Now want next placement.  Have 'nrem' rows still to place.
      do 900 irow = 1, n
         if (ifin(irow).eq.1) then
c           For all still-to-be-placed rows...
c           1. Place right at beg.:
               sim1 = 0.0
               do 500 j = 1, m
                  sim1 = sim1 + a(irow,j)*b(1,j)
  500          continue
c           2. Place right at end:
               if (nplace.gt.1) then
                  sim2 = 0.0
                  do 600 j = 1, m
                     sim2 = sim2 + a(irow,j)*b(nplace,j)
  600             continue
               endif
c           3. Place between k and k+1, where k = 1, ..., nplace-1:
               if (nplace.ge.2) then
                  do 800 k = 1, nplace-1
c                    Path length involves sim with k'th and with k+1'th rows
c                    in 'b'; i.e. b(k,j) and b(k+1,j), for all j.
c                    Sim is with a(irow,j), for all j.
                     sim3 = 0.0
                     do 700 j = 1, m
                        sim3 = sim3 + a(irow,j)*(b(k,j)+b(k+1,j))
  700                continue
                     if (sim3.gt.sim) then
                        sim = sim3
                        insrt = k
                        iplrow = irow
                     endif
  800             continue
               endif
c
c              Scale up 'sim1' and 'sim2' relative to 'sim', since former
c              are based on one link only
               sim1 = sim1*2.0
               sim2 = sim2*2.0
c              Use 'sim' and 'insrt' to store final info on row to place.
               if (sim1.gt.sim) then
                  sim = sim1
                  insrt = 0
                  iplrow = irow
               endif
c              .ge. in following, to force tied case to end
               if (sim2.ge.sim) then
                  sim = sim2
                  insrt = nplace+1
                  iplrow = irow
               endif
            endif
 900     continue
c
c      So now, we want to make placement in location 'insrt+1'
c      1. This happens to be right at beginning:
       if (insrt.eq.0) then
c         Shift right
          do 1100 l = nplace+1, 2, -1
             ib(l)  = ib(l-1)
             do 1000 j = 1, m
                b(l,j) = b(l-1,j)
 1000        continue
 1100     continue
          do 1200 j = 1, m
             b(1,j) = a(iplrow,j)
 1200     continue
          ifin(iplrow) = 0
          nplace       = nplace + 1
          nrem         = nrem - 1
          ib(1)        = iplrow
          goto 1900
       endif
c
c      2. Placement happens to be right at end of all current placements:
       if (insrt.eq.nplace+1) then
c         Insert after all current placements.
          do 1300 j = 1, m
             b(nplace+1,j) = a(iplrow,j)
 1300     continue
          ifin(iplrow) = 0
          nplace       = nplace + 1
          nrem         = nrem -1
          ib(nplace)   = iplrow
          goto 1900
       endif
c
c      3. If we get to here, new placement is somewhere in the middle.
c      Shift rows 'insrt+1' to 'nplace', in 'b', right.
       do 1500 l = nplace+1, insrt+2, -1
          ib(l) = ib(l-1)
          do 1400 j = 1, m
             b(l,j) = b(l-1,j)
 1400     continue
 1500  continue
       do 1600 j = 1, m
          b(insrt+1,j) = a(iplrow,j)
 1600  continue
       nplace       = nplace + 1
       nrem         = nrem -1
       ifin(iplrow) = 0
       ib(insrt+1)   = iplrow
       goto 1900
c
 1900   continue
       if (nrem.ge.1) goto 400
c
       return
       end
                 
c------------------------------------------------------------------------------
      subroutine cbea(n,m,a,jstart,b,jb,jfin)
c     Reoder cols. using BEA, bond energy algorithm.
c     See references at beg. of routine 'rbea'.
      dimension a(n,m), b(n,m), jb(m), jfin(m)
c------------------------------------------------------------------------------
c     a(n,m)         input matrix, cols. of which are to be permuted
c     jstart         1st col. to be placed
c     b(n,m)         permuted cols. to be stored in this array
c     jb(m)          integer list giving permutation carried out
c     jfin(m)        book-keeping vector: is col. still active, or has it been 
c                    placed (resp. = 1, = 0).
c------------------------------------------------------------------------------
c
c     Flags to indicate if col. already chosen
      do 200 j = 1, m
 200     jfin(j) = 1
c
c     'nplace' cols. placed (anticipating!); 'nrem' cols. still to be placed.
      nplace = 1
      nrem   = m-1
c
c     Place 1st col.
      do 300 i = 1, n
         b(i,1)     = a(i,jstart)
         jfin(jstart) = 0
         jb(nplace) = jstart
 300  continue
c
 400  continue
      sim   = -100000.0
      insrt = 0
c     'nplace' cols. have been placed.
c     Now want next placement.  Have 'nrem' cols. still to place.
      do 900 jcol = 1, m
         if (jfin(jcol).eq.1) then
c           For all still-to-be-placed cols...
c           1. Place right at beg.:
               sim1 = 0.0
               do 500 i = 1, n
                  sim1 = sim1 + a(i,jcol)*b(i,1)
  500          continue
c           2. Place right at end:
               if (nplace.gt.1) then
                  sim2 = 0.0
                  do 600 i = 1, n
                     sim2 = sim2 + a(i,jcol)*b(i,nplace)
  600             continue
               endif
c           3. Place between k and k+1, where k = 1, ..., nplace-1:
               if (nplace.ge.2) then
                  do 800 k = 1, nplace-1
c                    Path length involves sim with k'th and with k+1'th cols.
c                    in 'b'; i.e. b(i,k) and b(i,k+1), for all i.
c                    Sim is with a(i,jcol), for all i.
                     sim3 = 0.0
                     do 700 i = 1, n
                        sim3 = sim3 + a(i,jcol)*(b(i,k)+b(i,k+1))
  700                continue
                     if (sim3.gt.sim) then
                        sim = sim3
                        insrt = k
                        jplcol = jcol
                     endif
  800             continue
               endif
c
c              Scale up 'sim1' and 'sim2' rel. to 'sim' since former are
c              based on only one link
               sim1 = 2.0*sim1
               sim2 = 2.0*sim2
c              Use 'sim' and 'insrt' to store final info. on col. to be placed.
               if (sim1.gt.sim) then
                  sim = sim1
                  insrt = 0
                  jplcol = jcol
               endif
c              .ge. in following, to force tied case to end
               if (sim2.ge.sim) then
                  sim = sim2
                  insrt = nplace+1
                  jplcol = jcol
               endif
            endif
 900     continue
c
c      So now, we want to make placement in location 'insrt+1'.
c      1. This happens to be right at beginning.
       if (insrt.eq.0) then
c         Shift right
          do 1100 l = nplace+1, 2, -1
             jb(l) = jb(l-1)
             do 1000 i = 1, n
                b(i,l) = b(i,l-1)
 1000        continue
 1100     continue
          do 1200 i = 1, n
             b(i,1) = a(i,jplcol)
 1200     continue
          jfin(jplcol) = 0
          nplace       = nplace + 1
          nrem         = nrem - 1
          jb(1)   = jplcol
          goto 1900
       endif
c
c      2. Placement happens to be right at end of already placed cols.
       if (insrt.eq.nplace+1) then
          do 1300 i = 1, n
             b(i,nplace+1) = a(i,jplcol)
 1300     continue
          jfin(jplcol) = 0
          nplace       = nplace + 1
          nrem         = nrem -1
          jb(nplace)   = jplcol
          goto 1900
       endif
c
c      3. New col. placement happens to be somewhere in the middle.
c      Shift cols. 'insrt+1' to 'nplace', in 'b', right.
       do 1500 l = nplace+1, insrt+2, -1
          jb(l)  = jb(l-1)
          do 1400 i = 1, n
             b(i,l) = b(i,l-1)
 1400     continue
 1500  continue
       do 1600 i = 1, n
          b(i,insrt+1) = a(i,jplcol)
 1600  continue
       nplace       = nplace + 1
       nrem         = nrem -1
       jfin(jplcol) = 0
       jb(insrt+1)   = jplcol
       goto 1900
c
 1900   continue
       if (nrem.ge.1) goto 400
c
       return
       end
c--------------------------------------------------------------------------
       subroutine energy(n,m,b,ener)
       dimension b(n,m)
c      Det. "bond energy" of array b
c      I.e. product of each elt. with its 4 nearest neighbors,
c      summed over all elts.
       ener = 0.0
c      Corner elts.
       ener = ener + b(1,1)*(b(1,2)+b(2,1))
       ener = ener + b(1,m)*(b(1,m-1)+b(2,m))                 
       ener = ener + b(n,1)*(b(n-1,1)+b(n,2))
       ener = ener + b(n,m)*(b(n-1,m)+b(n,m-1))
c      Next non-corner border elts.
       do 100 j = 2, m-1
          ener = ener + b(1,j)*(b(1,j-1)+b(1,j+1)+b(2,j))
          ener = ener + b(n,j)*(b(n,j-1)+b(n,j+1)+b(n-1,j))
 100   continue
       do 200 i = 2, n-1
          ener = ener + b(i,1)*(b(i-1,1)+b(i+1,1)+b(i,2))
          ener = ener + b(i,m)*(b(i-1,m)+b(i+1,m)+b(i,m-1))
 200   continue
c      Finally, all non-border elts.
       do 400 i = 2, n-1
          do 300 j = 2, m-1
             ener = ener + b(i,j)*(b(i-1,j)+b(i+1,j)+b(i,j-1)+b(i,j+1))
 300      continue
 400   continue
c
       return
       end
       
