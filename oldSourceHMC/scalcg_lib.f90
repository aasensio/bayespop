module scalcg_mod
use forward_mod
contains
!             ======================================      
!             ******  SCALCG  ******  SCALCG  ******     
!             ======================================           
!
!-------------------------------------------------------------------!
!                       Subroutine SCALCG                                    
!
!             Scaled Conjugate Gradient Algorithms for 
!                   Unconstrained Optimization
!
!                         Neculai Andrei
!                Research Institute for Informatics
!           Center for Advanced Modeling and Optimization
!            8-10, Averescu Avenue, Bucharest 1, Romania
!                     E-mail: nandrei@ici.ro
!                        voice: 666.58.70                       
!                                                                   
!
! !-----------------------------------------------------------------!
! | SCALCG is a subroutine dedicated to compute the minimizer of a  |
! | differentiable function with a large number of variables.       |
! | SCALCG implements a nonlinear scaled conjugate gradient         |
! | algorithm. (Please see the paper: "Scaled Conjugate Gradient    |
! | Algorithms for Unconstrained Optimization")                     |
! | This subroutine is accompanied by subroutine "LineSearch" which |
! | implements the Wolfe line search. Both these subroutines belong |
! | to SCALCG package.                                              |
! |                                                                 |
! | The user must provide a subroutine to evaluate the function     |
! | value and its gradient in a point. The name of this subroutine  |
! | is EVALFG. On the web there is a file "func.for" where over 70  |
! | Fortran expression of some test functions are included.         |
! | Additionally, the algebraic expression of these functions is    |
! | given in file "testuo.pdf".                                     |
! |                                                                 |
! | There are some facilities for the user to specify:              |
! | 1) The termination criterion.                                   |
! | 2) Convergence tolerance for gradient.                          |
! | 3) Convergence tolerance for function value.                    |
! | 4) Maximum number of iterations in SCALCG subroutine.           |
! | 5) Maximum number of iterations in LineSearch subroutine.       |
! | 6) Level of printing.                                           |
! |                                                                 |
! !-----------------------------------------------------------------!
! |                                                                 |
! | The calling sequence of SCALCG is:                              |
! |                                                                 |
! |     call scalcg(n,x,epsg,epsf,delta,maxiter,stoptest,fx,gnorm,  |
! |                 iter,irstart,igamma,fgcnt,lscnt,tetas,tetaa,    |
! |                 iprint, nexp)                                   |
! |                                                                 |
! |                                                                 |
! |Input parameters:                                                |
! |=================                                                |
! |n          (integer) number of variables.                        |
! |x          (double)  starting guess, length n. On output         |
! |                     contains the solution.                      |
! |epsg       (double)  convergence tolerance for gradient.         |
! |epsf       (double)  convergence tolerance for function.         |
! |delta      (double)  parameter used in anticipative selection    |
! |                     the scaling factor of gradient. It has a    |
! |                     value comparable with the function's        |
! |                     value in the current point.                 |
! |maxiter    (integer) maximum number of iterations.               |
! |stoptest   (integer) parameter for selection of stopping         |
! |                     criterion.                                  |
! |                     If stoptest = 1, then consider the test:    |
! |                       if(ginf .le. epsg)                        |
! |                     If stoptest = 2, then consider the test:    |
! |                       if(ginf .le. dmax1(epsg, epsf*ginfz))     |
! |                     If stoptest = 3, then the test is:          |
! |                       if(gnorm .le. epsg)                       |
! |                     If stoptest = 4, then the test is:          |
! |                       if(gnorm .le. epsg*dmax1(1.d0, dabs(fx))) | 
! |                     If stoptest = 5, then the test is:          |
! |                       if(ginf .le. epsg/sqrt(n))                |
! |                     If stoptest = 6, then consider the test:    |
! |                       if(ginf .lt. epsg .OR.                    |
! |                          dabs(alfa*gtd) .le. epsf*dabs(fx))     |
! |                        where:                                   |
! |                        ginf  = infinite norm of gradient g(xk), |
! |                        ginfz = infinite norm of gradient g(x0), |
! |                        gnorm = norm 2 of gradient g(xk).        |
! |tetas      (logical) if tetas = .true. then the spectral formula |
! |                     is used by the algorithm.                   |
! |tetaa      (logical) if tetaa = .true. then the anticipative     |
! |                     formula is used by the algorithm.           |
! |iprint     (integer) parameter for printing the iterations.      |
! |                     if iprint .eq. 1 then print:                |
! |                        (iter, fgcnt, alpha, function, gnorm)    |
! |                     if iprint .ne. 1 then no print.             |
! |nexp       (integer) parameter specifying the number of the      |
! |                     problem considered in a train of experiments|
! |                                                                 |
! |                                                                 |
! |Output parameters:                                               |
! |==================                                               |
! |fx         (double)  function value in final (optimal) point.    |
! |gnorm      (double)  norm-2 of gradient in final point.          |
! |iter       (integer) number of iterations to get the final point.|
! |irstart    (integer) number of restart iterations.               |
! |igamma     (integer) number of gamma negative in anticipative    |
! |                     formula.                                    |
! |fgcnt      (integer) number of function evaluations.             |
! |lscnt      (integer) number of line searches.                    |
! !-----------------------------------------------------------------!
! |                                                                 |
! |                                                                 |
! |Calling subroutines:                                             |
! |====================                                             | 
! |Subroutine SCALCG is calling two subroutines:                    |
! |EVALFG       an user subroutine, and                             |
! |LINESEARCH   a package subroutine.                               |
! |                                                                 |   
! |The user must supply a subroutine with the function and its      |
! |gradient:                                                        |
! |  call evalfg(n,x,fx,grad,nexp)                                  |
! |where:                                                           |
! |  n    (integer)  number of variables.                           |
! |  x    (double)   the current iteration.                         |
! |  fx   (double)   function value in point x.                     |
! |  grad (double)   array with gradient of function in point x.    |
! |  nexp (integer)  parameter specifying the number of the         |
! |                  problem considered in a train of experiments.  |
! |                                                                 |
! |                                                                 |
! |The Wolfe line search is implemented in the subroutine           |
! |LINESEARCH, which belongs to the package.                        |
! |The calling sequence of this subroutine is as follows:           |            
! | call LineSearch (n,x,f,d,gtd,dnorm,alpha,xnew,fnew,gnew,        |
! |                  fgcnt,lscnt,lsflag, nexp)                      |
! |where:                                                           |
! |  n      (integer)  number of variables.                         |
! |  x      (double)   the current iteration.                       |
! |  f      (double)   function value in current point.             |
! |  d      (double)   array with search direction.                 |
! |  gtd    (double)   scalar: grad'*d.                             |   
! |  dnorm  (double)   2 norm of d.                                 |
! |  alpha  (double)   step length (given by the LineSearch).       |
! |  xnew   (double)   array with the new estimation of variables.  |
! |  fnew   (double)   function value in xnew.                      |
! |  gnew   (double)   array with gradient in xnew.                 |
! |  fgcnt  (integer)  number of function evaluations.              |
! |  lscnt  (integer)  number of line searches.                     |
! |  lsflag (integer)  parameter for abnormal Line Search           |
! |                    Termination. If the # of iterations in       |
! |                    LineSearch is greater than 20 then lsflag=1  |
! |  nexp   (integer)  parameter specifying the number of the       |
! |                    problem considered in a train of experiments |
! |                                                                 |
! |Subroutine LINESEARCH is the same as that used in CONMIN package |
! |by Shanno and Phua, and in SCG package by Birgin and Martinez.   |
! |                                                                 |
! !-----------------------------------------------------------------!
! |                                          Neculai Andrei, 2006   |
! !-----------------------------------------------------------------!
           
      subroutine scalcg(n,x,epsg,epsf,delta,maxiter,stoptest,fx,gnorm,&
                       iter,irstart,igamma,fgcnt,lscnt,tetas,tetaa, &
                       iprint,node)  
!
      parameter(ia=10000) 
      
      integer fgcnt,lscnt,lsflag
      integer stoptest, iprint
      integer n, node, iter, maxiter, irstart, igamma
      real*8  epsg, epsf, epsinf
      real*8 x(ia), xnew(ia), grad(ia), gradnew(ia), d(ia)
      real*8 s(ia),y(ia), ss(ia), ys(ia), v(ia),w(ia)
      real*8 fx, fxnew, gtg, g1tg, gtd, gnorm, dnorm, dnormpre
      real*8 alfa, theta, thetas, gtdr, dtd, eta, delta, alfan
      real*8 gts,gty,yts,yty,yks,yky,gtw,ytw, coef, sts,sty
      real*8 gamma, ginf, ginfz, gts1,yts1                
      logical tetas, tetaa

      
! Step 1. Initialization.
!
      epsinf = epsg/sqrt(float(n))
             
      iter=0   
      irstart=0           
      igamma=0
      
      call evalfg(n,x,fx,grad,nexp) 
      fgcnt=fgcnt+1

      gtg=0.d0
      ginfz=0.d0
      do i=1,n
        d(i) = -grad(i)
        gtg = gtg + grad(i)**2
        ginfz = dmax1(ginfz, dabs(grad(i)))
      end do 

      gtd = -gtg
      gnorm = sqrt(gtg)
      dnorm = gnorm
      gtdr = gtd
      dtd  = gtg

            
! * Step 2. Line search from initial point. Determine the steplength alfa.
! *         Initialy, alfa = 1/norm of gradient(x0).
! *         Compute: xnew = x + alfa*d, fnew=f(xnew), gradnew = gradient(xnew)

      if(gnorm .ne. 0.0d0) alfa=1.d0/gnorm  
        
      call linesearch(n,x,fx,d,gtd,dnorm,alfa, xnew,fxnew,gradnew,&
                     fgcnt,lscnt,lsflag, nexp)       


! *------------------------------------------------------ Write in file *.inf
          if(iprint .eq. 1) then
            if(mod(iter,10) .eq. 0) then
!               write(9, 330)
				  write(*, 330)
330           format(3x,'iter',3x,'fgcnt',10x,'alpha',18x,'function',&
                    16x,'gnorm')
            end if
            write(9,335) node, iter,fgcnt,alfa,fxnew,gnorm
				write(*,335) node, iter,fgcnt,alfa,fxnew,gnorm
335         format(2x,i2,2x,i5,3x,i5,3x,e20.13,3x,e20.13,3x,e20.13)
          end if                          
! *--------------------------------------------------------------------------


        
! * Step 3 Compute s and y. Prepare test for continuation 
! *                         ginf  = infinite norm of gradient
! *                         gtg   = gradnew'*gradnew (scalar product)
! *                         gnorm = norm 2 of gradient        
        gtg=0.d0
        ginf=0.d0     
      do i=1,n
        s(i) = xnew(i) - x(i)
        y(i) = gradnew(i) - grad(i) 
        gtg  = gtg + gradnew(i)*gradnew(i)
        ginf = max(ginf, dabs(gradnew(i)))
      end do  
     
      gnorm = sqrt(gtg)               
              
! * Step 4. Stop test for initial point. 

         if(gnorm .le. epsg*max(1.d0, abs(fxnew))) go to 999
      
! * Step 5. Continue iterations. 

         iter = iter + 1           

111   continue
! *------------------------------- Here a restart iteration,
! *                                i.e. the Powell restart criterion has been satisfied
! *
! * Step 6. Compute theta 
! *                                            |---------------------|
! *------------------------------------------- | 1) Spectral Formula |
! *                                            |---------------------|
      
      if(tetas) then
        sts = 0.d0
        sty = 0.d0
        do i=1,n
          sts = sts + s(i)*s(i)
          sty = sty + s(i)*y(i)
        end do
                        
        if(sty .le. 0.d0) then
          theta = 1.0d+10
        else
          theta = min(1.0d+10, max(1.0d-10, sts/sty))
        end if
      end if           
                          
! *                                            |-------------------------|
! *------------------------------------------- | 2) Anticipative Formula | 
! *                                            |-------------------------|
      
      if(tetaa) then
      
! c Here an alternative criterion for stopping the iterations:                                                           
! c          if( dabs(fxnew-fx)/(1.d0+dabs(fx)) .le. epsf) then
! c            istopfr=1
! c            go to 999
! c          end if  
                                                            
          gamma = 2.d0*(fxnew - fx - alfa*gtdr)/(dtd*alfa*alfa)   

        if(gamma .le. 0.d0) then       
          igamma = igamma+1
          eta   = (fx-fxnew+alfa*gtdr+delta)/gtdr  
          alfan = alfa - eta
          gamma = 2.d0*(fxnew - fx - alfan*gtdr)/(dtd*alfan*alfan)     
        end if 
        
        theta = 1.d0/gamma
      end if

! *------------------------------------------------ End theta computation
! *
! * Step 7. Update function value and norm of direction
                          
      fx=fxnew
      dnormpre = dnorm
        
                  
! * Step 8. Compute direction d.
! *                Here 4 scalar products for restart d(k+1) 
! *                Firstly update x and grad 
! *         Line search preparations for this new direction
! *                Compute gtd and dtd
! *
      gts=0.d0
      gty=0.d0
      yts=0.d0
      yty=0.d0                             
      do i=1,n
          x(i) = xnew(i)
          grad(i) = gradnew(i)
        gts = gts + grad(i)*s(i)
        gty = gty + grad(i)*y(i)
        yts = yts + y(i)*s(i)
        yty = yty + y(i)*y(i)
      end do
                 
      if(yts .ne. 0.d0) then        
        coef = 1.d0 + theta*yty/yts
        gtd=0.d0
        dtd=0.d0        
        do i=1,n
          d(i) = -theta * grad(i) &
                +theta * (gts/yts) * y(i)  &
                -( coef*gts/yts - theta*gty/yts ) * s(i)    
          gtd = gtd + grad(i)*d(i)    
          dtd = dtd + d(i)*d(i)          
        end do
      end if
      
      if(yts .eq. 0.d0) then
        gtd=0.d0
        dtd=0.d0
        do i=1,n
          d(i) = -theta * grad(i)
          gtd = gtd + grad(i)*d(i)    
          dtd = dtd + d(i)*d(i)
        end do
      end if      
            
      dnorm = sqrt(dtd)
      gtdr = gtd


! * Step 9. Line search.
! *         Initialize alfa and then with this value compute the step length
! *         Compute: xnew = x + alfa*d, fnew=f(xnew), gradnew = gradient(xnew)

      alfa = alfa * (dnormpre/dnorm)
      
      call linesearch(n,x,fx,d,gtd,dnorm,alfa, xnew,fxnew,gradnew,&
                     fgcnt,lscnt,lsflag, nexp)       

! *------------------------------------------------------ Write in file *.inf
          if(iprint .eq. 1) then
            if(mod(iter,10) .eq. 0) then
!               write(9, 330)
				  write(*, 330)
            end if
            write(9,335) node,iter,fgcnt,alfa,fxnew,gnorm
				write(*,335) node,iter,fgcnt,alfa,fxnew,gnorm
          end if                          
! *--------------------------------------------------------------------------
      

! * Step 10. Store theta.
                   
      thetas=theta
                         
! * Step 11. Compute s and y.  Store these vectors for normal computation.
! *          Prepare for stopping test.                  

        gtg  = 0.d0
        g1tg = 0.d0
        ginf = 0.d0
      do i=1,n       
        s(i) = xnew(i) - x(i)
        y(i) = gradnew(i) - grad(i)      
        ss(i) = s(i)
        ys(i) = y(i)        
        gtg = gtg + gradnew(i)*gradnew(i)
        g1tg = g1tg + gradnew(i)*grad(i)
        ginf = max(ginf, dabs(gradnew(i))) 
      end do  
      gnorm = sqrt(gtg)
                         
! * Step 12. Test of Stopping the iterations  ------------------------ (I) 
! *                                  
      if(stoptest .eq. 1) then
        if(ginf .le. epsg) go to 999 
      end if

      if(stoptest .eq. 2) then
        if(ginf .le. dmax1(epsg, epsf*ginfz)) go to 999
      end if    

      if(stoptest .eq. 3) then
        if(gnorm .le. epsg) go to 999
      end if    

      if(stoptest .eq. 4 ) then
        if(gnorm .le. epsg * dmax1(1.d0, dabs(fxnew))) go to 999
      end if    

      if(stoptest .eq. 5 ) then
        if(ginf .le. epsinf) go to 999
      end if    

      if(stoptest .eq. 6 ) then
        if(ginf .lt. epsg .OR. &
          dabs(alfa*gtd) .le. epsf*dabs(fxnew)) go to 999
      end if
! *-----------------------------------------------------------------------

! * Step 13. Continue iterations.

222   continue
 
      iter = iter + 1
      if(iter .gt. maxiter) then
! c        write(*,100)
! c100     format(4x,'Too many iterations')
        go to 1000
      end if                      
      

! * Step 14. Test for restart.             **** Powell criterion. ****             
! *                                        ===========================
! c        g1tg=0.d0
! c        do i=1,n
! c          g1tg = g1tg + gradnew(i)*grad(i)            
! c        end do
      
        if(dabs(g1tg) .ge. 0.2d0 * gtg) then 
          irstart = irstart + 1
          go to 111
        end if  


! * Step 15.       
! *------------------------------ Here Normal iteration.
! *                               The Powell restart criterion is not satisfied.
! *                               Compute the direction using the BFGS updating philosophy.              
      fx=fxnew           
      dnormpre=dnorm           


! * Step 16. Compute v, w, d
! *          Here 6 scalar products for v and w                    
! *               2 scalar products for d(k+1): gts1, yts1

      gts=0.d0
      gty=0.d0
      yts=0.d0
      yty=0.d0
      yks=0.d0
      yky=0.d0
        gts1=0.d0
        yts1=0.d0
      do i=1,n 
          x(i)    = xnew(i)
          grad(i) = gradnew(i)      
        gts = gts + grad(i)*ss(i)
        gty = gty + grad(i)*ys(i)
        yts = yts + ys(i)*ss(i)
        yty = yty + ys(i)*ys(i)
        yks = yks + y(i) *ss(i)
        yky = yky + y(i) *ys(i)
          gts1 = gts1 + grad(i)*s(i)
          yts1 = yts1 + y(i)*s(i)
      end do       
      
      coef = 1.d0 + thetas*yty/yts

        gtw=0.d0
        ytw=0.d0
      
      do i=1,n
        v(i) = thetas*grad(i) &
             -thetas*(gts/yts)*ys(i) &
             +(coef*gts/yts - thetas*gty/yts)*ss(i)
     
        w(i) = thetas*y(i)&
             -thetas*(yks/yts)*ys(i) &
             +(coef*yks/yts - thetas*yky/yts)*ss(i)

        gtw = gtw + grad(i)*w(i)
        ytw = ytw + y(i)   *w(i)
     
      end do
! *
! *            Here compute d(k+1).
! * Step 17.   Prepare for line search for this direction: gtd, dtd
! *                                 
      gtd=0.d0
      dtd=0.d0
      do i=1,n
        d(i) = -v(i)&
              +(gts1*w(i) + gtw*s(i))/yts1 &
              -(1.d0+ytw/yts1)*(gts1/yts1)*s(i)
     
        gtd = gtd + grad(i)*d(i)
        dtd = dtd + d(i)*d(i)
      end do

      dnorm = sqrt(dtd)
      gtdr =gtd


! * Step 18. Line search.
! *          Initialize alpha and then with this value compute the step length.
! *          Compute: xnew = x + alfa*d, fnew=f(xnew), gradnew = gradient(xnew)

      alfa = alfa * (dnormpre/dnorm)
      
      call linesearch(n,x,fx,d,gtd,dnorm,alfa, xnew,fxnew,gradnew,&
                     fgcnt,lscnt,lsflag, nexp)       


! *------------------------------------------------------ Write in file *.inf
          if(iprint .eq. 1) then
            if(mod(iter,10) .eq. 0) then
!               write(9, 330)
				  write(*, 330)
            end if
            write(9,335) node, iter,fgcnt,alfa,fxnew,gnorm
				write(*,335) node, iter,fgcnt,alfa,fxnew,gnorm
          end if                          
! *--------------------------------------------------------------------------
         
   
! * Step 19. Compute s and y.
! *          Prepare for test of continuation.                       

      gtg  = 0.d0
      g1tg = 0.d0
      ginf = 0.d0                                
      do i=1,n
          s(i) = xnew(i) - x(i)
          y(i) = gradnew(i) - grad(i)      
        gtg = gtg + gradnew(i)*gradnew(i) 
        g1tg = g1tg + gradnew(i)*grad(i)
        ginf =max(ginf, dabs(gradnew(i)))
      end do
      gnorm = sqrt(gtg)
                                

! * Step 20 Test of Stopping the iterations  --------------------- (II)

      if(stoptest .eq. 1) then
        if(ginf .le. epsg) go to 999 
      end if

      if(stoptest .eq. 2) then
        if(ginf .le. dmax1(epsg, epsf*ginfz)) go to 999
      end if    

      if(stoptest .eq. 3) then
        if(gnorm .le. epsg) go to 999
      end if    

      if(stoptest .eq. 4 ) then
        if(gnorm .le. epsg * dmax1(1.d0, dabs(fxnew))) go to 999
      end if    

      if(stoptest .eq. 5 ) then
        if(ginf .le. epsinf) go to 999
      end if    

      if(stoptest .eq. 6 ) then
        if(ginf .lt. epsg .OR. &
          dabs(alfa*gtd) .le. epsf*dabs(fxnew)) go to 999
      end if
! *-----------------------------------------------------------------------


! * Step 21. Continue iterations.

      go to 222      
      
! * Here criteria of STOP are satisfied.

999   continue       

! * Write down the solution, if the case.
                        
! *      write(*,802)(i,x(i),i=1,n)
! *802   format(10x,i5,4x,e20.13)       


1000  continue

      return
      end subroutine scalcg
! *-------------------------------------------------- End SCALCG            
      

    

! *             =================================  
! *             ******  WOLFE LINE SEARCH  ******       
! *             =================================      
      
! c------------------------------------------------------------------

      subroutine LineSearch (n,x,f,d,gtd,dnorm,alpha,xnew,fnew,gnew,&
                            fgcnt,lscnt,lsflag, nexp)

! C This is the one-dimensional line search used in CONMIN by Shanno and Phua.
! C Please see the subroutine LineSearch in SCG package by Birgin and Martinez.

! C     SCALAR ARGUMENTS
      integer n,fgcnt,lscnt,lsflag
      double precision f,gtd,dnorm,alpha,fnew

! C     ARRAY ARGUMENTS
      double precision x(n),d(n),xnew(n),gnew(n)

! C     LOCAL SCALARS
      integer i,lsiter
      double precision alphap,alphatemp,fp,dp,gtdnew,a,b

      lsflag = 0   

! * Maximum number of LineSearch is max$ls (now=20)
      
      max$ls=20
      
      alphap = 0.0d0
      fp     = f
      dp     = gtd

      do i = 1,n
          xnew(i) = x(i) + alpha * d(i)
      end do

! c*1
      call evalfg(n,xnew,fnew,gnew, nexp)
      fgcnt = fgcnt + 1

      gtdnew = 0.0d0
      do i = 1,n
          gtdnew = gtdnew + gnew(i) * d(i)
      end do

      lsiter = 0

! * Test whether the Wolfe line search conditions hve been met.

 10   if ( alpha * dnorm .gt. 1.0d-30 .and. lsiter .lt. max$ls .and.&
          .not. ( gtdnew .eq. 0.0d0 .and. fnew .lt. f ) .and.&
          ( ( fnew .gt. f + 1.0d-04 * alpha * gtd .or. &
          abs( gtdnew / gtd ) .gt. 0.9d0 ) .or. ( lsiter .eq. 0 .and.&
          abs( gtdnew / gtd ) .gt. 0.5d0 ) ) ) then

! * Thest whether the new point has a negative slope and a higher function
! * value than that corresponding to alpha=0. In this case, the search has 
! * passed through a local max and is heading for a distant local minimum.
! * Reduce alpha, compute the new corresponding point, its function value 
! * and gradient, as well as gtdnew.
! * Repeat this test until a good point has been found.

 20       if ( alpha * dnorm .gt. 1.0d-30 .and. fnew .gt. f .and.&
              gtdnew .lt. 0.0d0 ) then

              alpha = alpha / 3.0d0

              do i = 1,n
                  xnew(i) = x(i) + alpha * d(i)
              end do

! c*2
              call evalfg(n,xnew,fnew,gnew, nexp)
              fgcnt = fgcnt + 1

              gtdnew = 0.0d0
              do i = 1,n
                  gtdnew = gtdnew + gnew(i) * d(i)
              end do

              alphap = 0.0d0
              fp     = f
              dp     = gtd

              goto 20

          end if

! * Cubic interpolation to find a new trial point corresponding to alphatemp.

          a = dp + gtdnew - 3.0d0 * ( fp - fnew ) / ( alphap - alpha )
          b = a ** 2 - dp * gtdnew

          if ( b .gt. 0.0d0 ) then
              b = sqrt( b )
          else
              b = 0.0d0
          end if

          alphatemp = alpha - ( alpha - alphap ) * ( gtdnew + b - a ) /&
                     ( gtdnew - dp + 2.0d0 * b )

! * Test whether the line minimum has been bracketed.

          if ( gtdnew / dp .le. 0.0d0 ) then

! * Here the minimum has been bracketed.
! * Test whether the trial point lies sufficiently within the bracketed 
! * interval. If it does not, choose alphatemp as the midpoint of the
! * interval.

              if ( 0.99d0 * max( alpha, alphap ) .lt. alphatemp .or.&
                 alphatemp .lt. 1.01d0 * min( alpha, alphap ) ) then
                  alphatemp = ( alpha + alphap ) / 2.0d0
              end if

          else

! * Here the minimum has not been bracketed.
! * The trial point is too small, double the largest prior point.

              if ( gtdnew .lt. 0.0d0 .and. &
                 alphatemp .lt. 1.01d0 * max( alpha, alphap ) ) then
                  alphatemp = 2.0d0 * max( alpha, alphap )
              end if

! * The trial point is too large, halve the smallest prior point.

              if ( ( gtdnew .gt. 0.0d0 .and.&
                 alphatemp .gt. 0.99d0 * min( alpha, alphap ) ) .or.&
                 alphatemp .lt. 0.0d0 ) then
                  alphatemp = min( alpha, alphap ) / 2.0d0
              end if

          end if

! * Save and continue the search.

          alphap = alpha
          fp     = fnew
          dp     = gtdnew

          alpha = alphatemp

          do i = 1,n
              xnew(i) = x(i) + alpha * d(i)
          end do

! c*3
          call evalfg(n,xnew,fnew,gnew, nexp)
          fgcnt = fgcnt + 1

          gtdnew = 0.0d0
          do i = 1,n
              gtdnew = gtdnew + gnew(i) * d(i)
          end do

          lsiter = lsiter + 1

      goto 10                                                       
      
! c------------------------------ End if for stopping criteria
      end if

      if (lsiter .ge. max$ls) then
          lsflag = 1
      end if

      if (lsiter .ne. 0) then
          lscnt = lscnt + 1
      end if

      return

      end subroutine LineSearch
! *-------------------------------------------------- End of LineSearch
                                



! *             ==================================  
! *             ******  EXETIME Subroutine  ******       
! *             ==================================  
                                
! *----------------------------------------------------------------
! *  Date created       : May 30, 1995
! *  Date last modified : May 30, 1995
! *
! *  Subroutine for execution time computation.
! *
! *----------------------------------------------------------------
! *
	subroutine exetim(tih,tim,tis,tic, tfh,tfm,tfs,tfc)

	  integer*4 tih,tim,tis,tic
	  integer*4 tfh,tfm,tfs,tfc

	  integer*4 ti,tf
	  integer*4 ch,cm,cs
	  data ch,cm,cs/360000,6000,100/

	  ti=tih*ch+tim*cm+tis*cs+tic
	  tf=tfh*ch+tfm*cm+tfs*cs+tfc
	  tf=tf-ti
	  tfh=tf/ch
	  tf=tf-tfh*ch
	  tfm=tf/cm
	  tf=tf-tfm*cm
	  tfs=tf/cs
	  tfc=tf-tfs*cs

	  return
	end subroutine exetim
! *-------------------------------------------------- End of EXETIM
end module scalcg_mod