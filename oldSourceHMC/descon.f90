module descon_m
use like_m, only : slikelihoodNoNormalization

	real(kind=8) :: epsg, epsf, delta, fx, gnorm
	integer, parameter :: iisize=300, idsize=300, icsize=30
	integer :: iidata(iisize)
	real(kind=8) :: ddata(idsize)
	character(len=40) :: cdata(icsize)
	integer :: stoptest, iter, irstart, igamma, fgcnt, lscnt, nexp, iprint
	logical :: tetas, tetaa
	integer :: maxfg
	real(kind=8) :: dc, cc
	integer :: irs, ibeta, itheta
	logical :: angle, powell
	
contains

!-----------------------------------------------------------------------
! Solve the optimization problem using DESCON
!-----------------------------------------------------------------------
	subroutine maximumLikelihoodDescon(sdim, x)
	integer :: sdim
	real(kind=8) :: x(sdim)

		epsg = 1.d-7
		epsf = 1.d-10
		delta = 1.d0
		maxiter = 200
		stoptest = 1
		iprint = 1

		tetas=.true.
		tetaa=.true.
		iprint = 1
		fgcnt = 0
		lscnt = 0
		
		dc = 7.d0/8.d0         
		cc = 0.5d0		
		angle = .false.
		powell = .TRUE.
		maxfg = 150

! Optimize the merit function using DESCON			
		call descon(sdim,x,epsg,maxiter,maxfg,fx,gnorm,stoptest,dc,&
			cc,iter,irs,fgcnt,lscnt,angle, powell,     &
			ibeta, itheta, nexp)
	
	end subroutine maximumLikelihoodDescon


!-----------------------------------------------------------------------
! Evaluate forward problem
!-----------------------------------------------------------------------
	subroutine evalfg(n,x,f,g,nexp)

	integer :: n, nexp
   real(kind=8) :: x(n), f, g(n)
		call slikelihoodNoNormalization(x,f,g,.TRUE.)
		f = -f
		g = -g
	
	end subroutine evalfg

!            *** Conjugate Gradient Algorithms Project *** 
!           ===============================================
!
!                                         
!
!
!-------------------------------------------------------------------
!                           Subroutine DESCON                                
!                         =====================
!                                   
!
!                             Neculai Andrei
!                  Research Institute for Informatics
!             Center for Advanced Modeling and Optimization
!              8-10, Averescu Avenue, Bucharest 1, Romania
!                        E-mail: nandrei@ici.ro
!                           voice: 666.58.70
!                                and
!                     Academy of Romanian Scientists 
!              Science and Information Technology Section              
!            54, Splaiul Independentei, Bucharest 5, Romania
!
!
!
! /-----------------------------------------------------------------\
! | DESCON is a subroutine dedicated to compute the minimizer of    |
! | a differentiable function with a large number of variables.     |
! |                                                                 |
! | This subroutine is accompanied by subroutine "LineSearch" which |
! | implements the Wolfe line search. Both these subroutines belong |
! | to DESCON package.                                              |
! |                                                                 |
! | The user must provide a subroutine to evaluate the function     |
! | value and its gradient in a point. The name of this subroutine  |
! | is EVALFG.                                                      |
! | The algebraic expression of the functions considered in this    |
! | program, as well as their Fortran code are presented into the   |
! | paper:                                                          |
! | N. Andrei, "An unconstrained optimization test functions        |
! | collection", Advanced Modeling and Optimization, vol.10, No.1,  |
! | (2008) pp.147-161.                                              | 
! | http://www.ici.ro/camo/journal/v10n1.htm                        |
! |                                                                 |
! | There are some facilities for the user to specify:              |
! | 1) The termination criterion.                                   |
! | 2) Convergence tolerance for gradient.                          |
! | 3) Convergence tolerance for function value.                    |
! | 4) Maximum number of iterations in LineSearch subroutine.       |
! | 5) The parameter w from the descent condition.                  |
! | 6) The parameter v from the conjugacy condition.                |
! |                                                                 |
! |-----------------------------------------------------------------|      
! |                                                                 |
! | The calling sequence of DESCON is:                              |
! |                                                                 |
! |    subroutine descon(n,x,epsg,maxiter,maxfg,f,gnorm,stoptest,dc,|
! |                      cc,iter,irs,fgcnt,lscnt,angle, powell,     |
! |                      ibeta, itheta, nexp)                       |
! |                                                                 |
! |Input parameters:                                                |
! |=================                                                |
! |n          (integer) number of variables.                        |
! |x          (double)  starting guess, length n. On output         |
! |                     contains the solution.                      |
! |epsg       (double)  convergence tolerance for gradient.         |
! |maxiter    (integer) maximum number of iterations.               |
! |maxfg      (integer) maxumum number of function and its gradient |
! |                     evaluations.                                |             
! |stoptest = option parameter for selection of                     | 
! |            stopping criterion:                                  |
! |            if stoptest = 1 then consider the following test:    | 
! |               if(ginf .le. epsg)                                |
! |            if stoptest = 2 then consider the following test:    |
! |               if(gnorm .le. epsg)                               |
! |               where:                                            |
! |               ginf  = infinite norm of gradient g(xk),          | 
! |               gnorm = norm-2 of gradient g(xk).                 | 
! |dc        (double)  parameter for sufficient descent condition   |
! |cc        (double)  parameter for conjugacy condition            |
! |angle      (logical) parameter specifying the angle criterion of |
! |                     restart.                                    |
! |powell     (logical) parameter specifying the Powell criterion of|
! |                     restart.                                    |   
! |nexp       (integer) parameter specifying the number of the      |
! |                     problem considered in a train of experiments|
! |                                                                 |
! |                                                                 |
! |Output parameters:                                               |
! |==================                                               |
! |f          (double)  function value in final (optimal) point.    |
! |gnorm      (double)  norm-2 of gradient in final point.          |
! |iter       (integer) number of iterations to get the final point.|
! |irs        (integer) number of restart iterations.               |
! |fgcnt      (integer) number of function evaluations.             |
! |lscnt      (integer) number of line searches.                    |
! |ibeta      (integer) number of iterations in which beta is zero  |
! |itheta     (integer) number of iterations in which theta is one  |
! |-----------------------------------------------------------------|
! |                                                                 |
! |                                                                 |
! |Calling subroutines:                                             |
! |====================                                             | 
! |Subroutine DESCON is calling two subroutines:                    |
! |EVALFG     an user subroutine (function and gradient),           |
! |LINESEARCH a package subroutine.                                 |
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
! | call LineSearch (n,x,f,d,gtd,dnorm,alpha,xnew,fnew,gnew,sigma,  |
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
! |  sigma  (double)   parameter sigma in the second Wolfe line     |
! |                    search condition. (input parameter)          |
! |  fgcnt  (integer)  number of function evaluations.              |
! |  lscnt  (integer)  number of line searches.                     |
! |  lsflag (integer)  parameter for abnormal Line Search           |
! |                    Termination. If the # of iterations in       |
! |                    LineSearch is greater than 20 then lsflag=1  |
! |  nexp   (integer)  parameter specifying the number of the       |
! |                    problem considered in a train of experiments |
! |                                                                 |
! |Subroutine LINESEARCH is the same as that used in CONMIN package |
! |by Shanno and Phua, and in SCALCG package by Andrei.             |
! |                                                                 |
! |-----------------------------------------------------------------|
! |                                           Neculai Andrei, 2011  |
! \-----------------------------------------------------------------/
!           
!
!
!********************************************************************

      subroutine descon(n,x,epsg,maxiter,maxfg,f,gnorm,stoptest,dc,cc,&
                       iter,irs,fgcnt,lscnt,angle, powell, &
                       ibeta, itheta, nexp)

      integer :: ia, n5, n6, ibeta, itheta
      parameter(ia=1000000)

!     SCALAR ARGUMENTS
      integer n,iter,irs,fgcnt,lscnt,maxiter,maxfg
      integer stoptest, nexp
      double precision epsg, f,gnorm     
      logical angle, powell 

!     ARRAY ARGUMENTS
      double precision x(n)

!     LOCAL SCALARS
      integer i,lsflag
      double precision fnew,alpha,gtg, beta, sts,&
            dnorm,dnormprev, ginf,&
            gtd, gtgp, stg, acc,bdc,  &
            yts, ytg, theta, delta, dc,cc, epsm, bn,sigma
      
!     LOCAL ARRAYS
      double precision xnew(ia),g(ia),gnew(ia),d(ia),&
             y(ia), s(ia)
      common /acca/epsm
            
! Initialization
      
      n5 = mod(n,5)
      n6 = n5 + 1
      
      iter   = 0
      irs    = 0
      fgcnt  = 0
      lscnt  = 0         
      ibeta  = 0
      itheta = 0
      
      call evalfg(n,x,f,g, nexp)
      fgcnt = fgcnt + 1

        gtg = 0.0d0
        do i = 1,n5
          d(i) = - g(i)
          gtg  = gtg + g(i) ** 2
        end do    
        do i = n6,n,5
          d(i)   = -g(i)
          d(i+1) = -g(i+1)
          d(i+2) = -g(i+2)
          d(i+3) = -g(i+3)
          d(i+4) = -g(i+4)
          gtg  = gtg + g(i)**2 + g(i+1)**2 + g(i+2)**2 +&
                                g(i+3)**2 + g(i+4)**2   
        end do
      gnorm = sqrt( gtg )

      gtd   = -gtg
      dnorm = gnorm                   
    

      if ( gnorm .gt. 0.0d0 ) then
          alpha = 1.0d0 / dnorm
      end if 

! Initial value of parameter sigma. 
! At the next iteration it is computed as in (3.11). Please see 
! ICI Technical Report, November 9, 2011.     

      sigma = 0.8d0       
    

    
! --------------------------------   Main loop   --------------------    
!====================================================================

110    continue
    
      
!------------------------------------  STOP test section
!                                      =================

      if(iter .eq. 0) go to 91

        if(stoptest .eq. 1) then
          ginf=dabs(g(1))
          do i=2,n5
            if(dabs(g(i))   .gt. ginf) ginf = dabs(g(i))
          end do   
          do i=n6,n,5
            if(dabs(g(i))   .gt. ginf) ginf = dabs(g(i))
            if(dabs(g(i+1)) .gt. ginf) ginf = dabs(g(i+1))
            if(dabs(g(i+2)) .gt. ginf) ginf = dabs(g(i+2))
            if(dabs(g(i+3)) .gt. ginf) ginf = dabs(g(i+3))
            if(dabs(g(i+4)) .gt. ginf) ginf = dabs(g(i+4))
          end do  
          if(ginf .le. epsg) go to 999    
        end if
!
        if(stoptest .eq. 2) then
          if(gnorm .le. epsg) go to 999
        end if 
!      
91    continue

      

!---------------------------------- Increment iteration section
!                                   ===========================

          iter = iter + 1  
          if(iter .gt. maxiter) go to 999         
          



!---------------------------------- Line search section
!                                   ===================                                           
!
! Determine the step length ALPHA and the new point XNEW, as well as
! the function value in xnew, FNEW, and the gradient in xnew, GNEW.


          call LineSearch(n,x,f,d,gtd,dnorm,alpha,xnew,fnew,gnew,&
                         sigma,fgcnt,lscnt,lsflag, n5,n6,nexp)

          if(fgcnt .gt. maxfg) go to 999

      alpha = alpha



!
!---------------------------------- Acceleration section
!                                   ====================
!
! (We use an/bn only if b is different from zero)
! (an=gtd)
      
      bn=0.d0
      do i = 1,n5       
        bn = bn + (g(i)-gnew(i))*d(i)
      end do   
      do i = n6,n,5
        bn= bn + (g(i)-gnew(i))*d(i)+&
                (g(i+1)-gnew(i+1))*d(i+1)+&
                (g(i+2)-gnew(i+2))*d(i+2)+&
                (g(i+3)-gnew(i+3))*d(i+3)+&
                (g(i+4)-gnew(i+4))*d(i+4)   
      end do
!
      if(dabs(bn) .gt. epsm) then
          do i=1,n5
            xnew(i)   = x(i)   + (gtd/bn)*alpha*d(i)
          end do   
          do i=n6,n,5
            xnew(i)   = x(i)   + (gtd/bn)*alpha*d(i)
            xnew(i+1) = x(i+1) + (gtd/bn)*alpha*d(i+1)
            xnew(i+2) = x(i+2) + (gtd/bn)*alpha*d(i+2)
            xnew(i+3) = x(i+3) + (gtd/bn)*alpha*d(i+3)
            xnew(i+4) = x(i+4) + (gtd/bn)*alpha*d(i+4)
          end do
          call evalfg(n,xnew,fnew,gnew, nexp) 
          fgcnt = fgcnt + 1
          if(fgcnt .gt. maxfg) go to 999
      end if


!
!
!---------------------------------- Prepare some scalar products
!                                   ============================
!              
         gtg = 0.d0
         gtgp= 0.d0          
         ytg = 0.d0 
         stg = 0.d0
         yts = 0.d0   
         sts = 0.d0
       do i = 1,n5    
         s(i) = xnew(i) - x(i)
         y(i) = gnew(i) - g(i)
         ytg = ytg + gnew(i)*y(i)
         stg = stg + gnew(i)*s(i)
         yts = yts + y(i)*s(i)
         gtgp = gtgp + gnew(i)*g(i)
           x(i) = xnew(i)
           g(i) = gnew(i)
         gtg  = gtg + g(i) * g(i)             
         sts  = sts + s(i) * s(i)
       end do  

       do i = n6,n,5
         s(i)   = xnew(i)   - x(i) 
         s(i+1) = xnew(i+1) - x(i+1)
         s(i+2) = xnew(i+2) - x(i+2)
         s(i+3) = xnew(i+3) - x(i+3)
         s(i+4) = xnew(i+4) - x(i+4)  
         y(i)   = gnew(i)   - g(i)
         y(i+1) = gnew(i+1) - g(i+1)
         y(i+2) = gnew(i+2) - g(i+2)
         y(i+3) = gnew(i+3) - g(i+3)
         y(i+4) = gnew(i+4) - g(i+4)
         ytg = ytg + gnew(i)*y(i)+gnew(i+1)*y(i+1)+gnew(i+2)*y(i+2)+&
                    gnew(i+3)*y(i+3)+gnew(i+4)*y(i+4)    
         stg = stg + gnew(i)*s(i)+gnew(i+1)*s(i+1)+gnew(i+2)*s(i+2)+&
                    gnew(i+3)*s(i+3)+gnew(i+4)*s(i+4)    
         yts = yts + y(i)*s(i)+y(i+1)*s(i+1)+y(i+2)*s(i+2)+&
                              y(i+3)*s(i+3)+y(i+4)*s(i+4)
         gtgp = gtgp + gnew(i)*g(i)+gnew(i+1)*g(i+1)+gnew(i+2)*g(i+2)+&
                     gnew(i+3)*g(i+3)+gnew(i+4)*g(i+4)
         sts = sts + s(i)*s(i)+s(i+1)*s(i+1)+s(i+2)*s(i+2)+&
                              s(i+3)*s(i+3)+s(i+4)*s(i+4)

           x(i)   = xnew(i)
           x(i+1) = xnew(i+1)
           x(i+2) = xnew(i+2)
           x(i+3) = xnew(i+3)
           x(i+4) = xnew(i+4)  
           g(i)   = gnew(i)
           g(i+1) = gnew(i+1)
           g(i+2) = gnew(i+2)
           g(i+3) = gnew(i+3)        
           g(i+4) = gnew(i+4)                                   
           gtg  = gtg + g(i)*g(i)+g(i+1)*g(i+1)+g(i+2)*g(i+2)+&
                       g(i+3)*g(i+3)+g(i+4)*g(i+4)      
       end do
!    
          gnorm= sqrt( gtg )                  
!
          f = fnew
          dnormprev = dnorm  
                
!                         
!--------------------------------- Sigma computation
!                                  =================
!
!
          sigma = gtg/(dabs(ytg)+gtg)
!  
!
!
! -------------------------------- Delta, Theta and Beta computation  
!                                  =================================
!
!1  Delta computation. (Delta bar)
!
          delta = stg*ytg-gtg*yts

!   acc (conjugacy condition)   (ak)
!   bdc (descent condition)     (bk)
                            
          acc = cc*stg + ytg
          bdc = dc*yts*gtg + ytg*stg
               
!
!2  Beta and theta computation                  
!   If |delta| is greater than epsilon machine, then compute
!   beta and theta by means of the algorithm.        
!   Otherwise, set theta=1 and beta = 0.
!
!   Now theta computation:

         if(dabs(delta) .gt. epsm) then 
          theta = (acc/ytg)*(1.d0+yts*gtg/delta)-bdc/delta
         else
          theta = 1.d0 
          itheta=itheta+1
         end if  
!
!
! Variants for beta computation
! =============================
!
! For beta computation we consider 3 possibilities:    
!
! 1) beta as in algorithm DESCON,
! 2) beta as suggested by Dai and Liao+ (i.e. similar to PRP+),
! 3) beta truncated (see Hager and Zhang).             
!
!
! Therefore, the algorithm can be implemented in 4 variants as follows:    
!
! a) (2.1),(2.2) and (2.12) for theta and (2.13) for beta computation.
!    
! b) (2.1),(2.2) and (2.12) for theta and (7.29) for beta computation.
!
! c) (2.1),(2.2) and (2.12) for theta and (2.13) for beta comptation 
!    with truncation.
!
! d) (2.1) (2.2) and (2.12) for theta and (7.29) for beta computation
!    with truncation.    
!
! Please see ICI Technical Report, November 9, 2011.
!
! Numerical experiments proved that all these variants have similar 
! performances.
! They do not improve significantly the performances of the algorithm.
! However, variant a) is slightly better.
!
!         
      if(yts .ne. 0.d0 .and. dabs(delta) .gt. epsm) then
       beta = ytg/yts-(bdc*ytg)/(yts*delta)+acc*gtg/delta
!c       beta = dmax1(ytg/yts,0.d0) -
!c     *        (bdc*ytg)/(yts*delta) +
!c     *        (acc*gtg)/delta
!c       beta = dmax1(beta,-1.d0/(dnorm*dmin1(0.1d0,gnorm)))
      else
        beta  = 0.d0              
        ibeta=ibeta+1
      end if  
!
!                        
! --------------------------------- Direction computation
!                                   =====================
!
         dnorm = 0.0d0
         gtd   = 0.0d0
         do i = 1,n5
           d(i)  = -theta*g(i) + beta*s(i) 
           dnorm =   dnorm + d(i) ** 2
           gtd   =   gtd + g(i) * d(i)  
         end do
         
         do i = n6,n,5  
           d(i)    =-theta * g(i)   + beta * s(i)     
           d(i+1)  =-theta * g(i+1) + beta * s(i+1) 
           d(i+2)  =-theta * g(i+2) + beta * s(i+2)   
           d(i+3)  =-theta * g(i+3) + beta * s(i+3)   
           d(i+4)  =-theta * g(i+4) + beta * s(i+4)   

           dnorm   =   dnorm + d(i)**2+d(i+1)**2+d(i+2)**2+&
                                      d(i+3)**2+d(i+4)**2
           gtd     =   gtd + g(i)*d(i)+g(i+1)*d(i+1)+g(i+2)*d(i+2)+&
                                      g(i+3)*d(i+3)+g(i+4)*d(i+4) 
         end do

         dnorm = sqrt( dnorm )   
    
!       
!------------------------------------ END direction computation
!             
!   RESTART CRITERIA
!   ================
!
! Angle Restart Test
!
        if(angle) then
          if ( gtd .gt. -1.0d-03 * gnorm * dnorm ) then
              irs = irs + 1
              do i = 1,n
                d(i) = -g(i)
              end do
              dnorm =  gnorm
              gtd   = -gtg  
          end if        
        end if

!
! Beale-Powell restart test
!             
        if(powell) then
          if(dabs(gtgp) .gt. 0.2d0*dabs(gtg)) then
              irs = irs + 1
              do i = 1,n5
                d(i) = -g(i) 
              end do     
              do i = n6,n,5
                d(i)   = -g(i)     
                d(i+1) = -g(i+1)   
                d(i+2) = -g(i+2)   
                d(i+3) = -g(i+3)   
                d(i+4) = -g(i+4)            
              end do              
              dnorm =   gnorm
              gtd   =  -gtg  
          end if            
        end if  
!------------------------------------------ Prepare first trial
!                                           of steplength          
!                                           ===================
          
        if(dnorm .ne. 0.d0) then
          alpha = alpha * dnormprev / dnorm
        else
          alpha = 1.d0
        end if    
!
      go to 110


!------------------------------------------ End of main loop
!                                           ================

999   continue
      
!     
      return 

      end subroutine descon

!-------------------------------------------- END DESCON subroutine





!******************************************************************

      subroutine LineSearch (n,x,f,d,gtd,dnorm,alpha,xnew,fnew,gnew,&
                            sigma,fgcnt,lscnt,lsflag, n5,n6,nexp)

                            
      integer :: n5, n6, max$ls
!     This is the one-dimensional line search used in CONMIN

!     SCALAR ARGUMENTS
      integer n,fgcnt,lscnt,lsflag, nexp
      double precision f,gtd,dnorm,alpha,fnew

!     ARRAY ARGUMENTS
      double precision x(n),d(n),xnew(n),gnew(n)

!     LOCAL SCALARS
      integer i,lsiter
      double precision alphap,alphatemp,fp,dp,gtdnew,a,b,sigma

      lsflag = 0                                              
      
! Maximum number of LineSearch is max$ls (now=6)     

      max$ls=8          
      
      alphap = 0.0d0
      fp     = f
      dp     = gtd

      do i = 1,n5
        xnew(i) = x(i) + alpha * d(i)
      end do     
      do i = n6,n,5
        xnew(i)   = x(i)   + alpha * d(i)
        xnew(i+1) = x(i+1) + alpha * d(i+1)
        xnew(i+2) = x(i+2) + alpha * d(i+2)
        xnew(i+3) = x(i+3) + alpha * d(i+3)
        xnew(i+4) = x(i+4) + alpha * d(i+4)
      end do    
!1
         call evalfg(n,xnew,fnew,gnew, nexp)
         fgcnt = fgcnt + 1

      gtdnew = 0.0d0
      do i = 1,n5
       gtdnew = gtdnew + gnew(i) * d(i)
      end do     
      do i = n6,n,5
       gtdnew = gtdnew + gnew(i)*d(i)+gnew(i+1)*d(i+1)+gnew(i+2)*d(i+2)+&
                        gnew(i+3)*d(i+3)+gnew(i+4)*d(i+4)
      end do
          
      lsiter = 0                                          

 10   if ( alpha * dnorm .gt. 1.0d-30 .and. lsiter .lt. max$ls .and.&
        .not. ( gtdnew .eq. 0.0d0 .and. fnew .lt. f ) .and.&
        ( ( fnew .gt. f + 1.0d-04 * alpha * gtd .or. &
        dabs( gtdnew / gtd ) .gt. sigma ) .or. ( lsiter .eq. 0 .and.&
        dabs( gtdnew / gtd ) .gt. 0.50d0 ) ) ) then

 20       if ( alpha * dnorm .gt. 1.0d-30 .and. fnew .gt. f .and.&
              gtdnew .lt. 0.0d0 ) then

              alpha = alpha / 3.0d0

              do i = 1,n5
                xnew(i) = x(i) + alpha * d(i)
              end do     
              do i = n6,n,5                
                xnew(i)   = x(i)   + alpha * d(i)
                xnew(i+1) = x(i+1) + alpha * d(i+1)
                xnew(i+2) = x(i+2) + alpha * d(i+2)
                xnew(i+3) = x(i+3) + alpha * d(i+3)
                xnew(i+4) = x(i+4) + alpha * d(i+4)
              end do  
!2
                call evalfg(n,xnew,fnew,gnew, nexp)
                fgcnt = fgcnt + 1

              gtdnew = 0.0d0
              do i = 1,n5
                gtdnew = gtdnew + gnew(i) * d(i)
              end do     
              do i = n6,n,5              
                gtdnew = gtdnew + gnew(i)*d(i)+gnew(i+1)*d(i+1)+&
                                 gnew(i+2)*d(i+2)+gnew(i+3)*d(i+3)+&
                                 gnew(i+4)*d(i+4)     
              end do

              alphap = 0.0d0
              fp     = f
              dp     = gtd

              goto 20

          end if
                 
          a = dp + gtdnew - 3.0d0 * ( fp - fnew ) / ( alphap - alpha )
          b = a ** 2 - dp * gtdnew

          if ( b .gt. 0.0d0 ) then
              b = sqrt( b )
          else
              b = 0.0d0
          end if

          alphatemp = alpha - ( alpha - alphap ) * ( gtdnew + b - a ) /&
                     ( gtdnew - dp + 2.0d0 * b )

          if ( gtdnew / dp .le. 0.0d0 ) then

              if ( 0.99d0 * max( alpha, alphap ) .lt. alphatemp .or.&
                 alphatemp .lt. 1.01d0 * min( alpha, alphap ) ) then
                  alphatemp = ( alpha + alphap ) / 2.0d0
              end if

          else

              if ( gtdnew .lt. 0.0d0 .and. &
                 alphatemp .lt. 1.01d0 * max( alpha, alphap ) ) then
                  alphatemp = 2.0d0 * max( alpha, alphap )
              end if

              if ( ( gtdnew .gt. 0.0d0 .and.&
                 alphatemp .gt. 0.99d0 * min( alpha, alphap ) ) .or.&
                 alphatemp .lt. 0.0d0 ) then
                  alphatemp = min( alpha, alphap ) / 2.0d0
              end if

          end if

          alphap = alpha
          fp     = fnew
          dp     = gtdnew

          alpha = alphatemp

          do i = 1,n5
            xnew(i) = x(i) + alpha * d(i)
          end do                           
          do i = n6,n,5
            xnew(i)   = x(i)   + alpha * d(i)
            xnew(i+1) = x(i+1) + alpha * d(i+1)
            xnew(i+2) = x(i+2) + alpha * d(i+2)
            xnew(i+3) = x(i+3) + alpha * d(i+3)
            xnew(i+4) = x(i+4) + alpha * d(i+4)
          end do   
!3
            call evalfg(n,xnew,fnew,gnew, nexp)
            fgcnt = fgcnt + 1

          gtdnew = 0.0d0
          do i = 1,n5
            gtdnew = gtdnew + gnew(i) * d(i)
          end do                              
          do i = n6,n,5
            gtdnew = gtdnew + gnew(i)*d(i)+gnew(i+1)*d(i+1)+&
                             gnew(i+2)*d(i+2)+gnew(i+3)*d(i+3)+     &
                             gnew(i+4)*d(i+4)
          end do

          lsiter = lsiter + 1

        goto 10

      end if

      if ( lsiter .ge. max$ls ) then
          lsflag = 1
      end if

      if ( lsiter .ne. 0 ) then
          lscnt = lscnt + 1
      end if

      return

      end subroutine LineSearch
!---------------------------------- End LineSearch subroutine

                                           
                                           

                                
                                
!-----------------------------------------------------------
!  Date created       : May 30, 1995
!  Date last modified : May 30, 1995
!
!  Subroutine for execution time computation.
!
!-----------------------------------------------------------
!
	subroutine exetim(tih,tim,tis,tic, tfh,tfm,tfs,tfc)
!
	  integer*4 tih,tim,tis,tic
	  integer*4 tfh,tfm,tfs,tfc
!
	  integer*4 ti,tf
	  integer*4 ch,cm,cs
	  data ch,cm,cs/360000,6000,100/
!
	  ti=tih*ch+tim*cm+tis*cs+tic
	  tf=tfh*ch+tfm*cm+tfs*cs+tfc
	  tf=tf-ti
	  tfh=tf/ch
	  tf=tf-tfh*ch
	  tfm=tf/cm
	  tf=tf-tfm*cm
	  tfs=tf/cs
	  tfc=tf-tfs*cs
!
	  return
	end subroutine exetim
	
end module descon_m
!---------------------------------------------- End of EXETIM

 
!---------------------------------- EVALFG subroutine - last line

!                                                  Neculai Andrei
!                                        Last line DESCON package
!================================================================
! Last Line  
