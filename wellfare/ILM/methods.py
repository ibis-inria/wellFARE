"""

This module implements the general functions for the linear inversion
methods, which are used in all the 1-step and 2-step methods

"""


import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fmin
from ..curves import Curve
    
import scipy.linalg as linalg
import numpy.linalg as linalg

try:
  import cvxopt
  cvxopt.solvers.options["show_progress"] = False
  cvm =  cvxopt.base.matrix
  CVXOPT_FOUND = True
except ImportError:
  CVXOPT_FOUND = False



# alphas (= regularization parameters) checked during cross-validation.
# They can be provided to the functions, by default we provide a large range of
# alphas so that it should work with the inference of growth rate, promoter
# activity, and protein concentration alike. 

DEFAULT_ALPHAS = np.logspace(-10, 10, 1000)


def make_iL(Nu, Nic, eps=.0001):
    """ Returns the inverse of a derivation matrix.
    
    Parameters
    -----------

    Nu
      Number of values of the signal u

    ic
      Number of other parameters (initial conditions) in
      the model

    eps_L
      Negligible factor.


    """

    ar = np.arange(Nu)
    iLu = 1.0*((ar - np.array([ar]).T)<=0)
    iLu[:,0] = (1.0/eps)

    if Nic==0:
      return iLu
    
    iLic = (1.0/eps)*np.identity(Nic)
    void = np.zeros((Nu, Nic))
    
    return np.vstack([np.hstack([iLic, void.T]),
                      np.hstack([void, iLu])])


def make_L(Nu, Nic, eps_L=.001):
    """ Returns a derivation matrix.

    Parameters
    -----------

    Nu
      Number of values of the signal u

    ic
      Number of other parameters (initial conditions) in
      the model

    eps_L
      Negligible factor.


    """


    L = np.identity(Nu+Nic)
    for i in range(Nic+1):
        L[i,i] = eps_L
    for i in range(Nic, Nu+Nic-1):
        L[i+1,i] = -1
    return L

#@profile
def GCV(y, H, alphas , Qv=None, optimize=False, send_state=None):
    """ Generalized Cross-Validation.
    
    Finds the right regularization parameter alpha
    for the following L2-Regularization problem:
    
      find x minimizing:
          \| Hx - y \|^2 + \alpha^2 \| x \|^2
    
    
    
    Parameters
    -----------
    
    y
      Vector of (noisy) observations.
    
    alphas
      Values of alpha to be tried.
    
    H
      The observation matrix of the problem, such that
      Hx = y, where y is a vector of observations, x is
      the vector to be estimated.
    
    Qv
      Optional, H can be provided instead.
      Tuple Q,v,Q2 with Q matrix and v vector of eigenvalues
      such that H = Q diag(v) Q.T, and Q2 = Q**2 (not dot)
      Note that Qv can be obtained from H with
    
      >>> Qv = np.linalg.eigh( H.dot(H.T) )
      >>> Q2 = Q**2
    
    
    Returns
    --------
    
    best_alpha, best_x, alphas_scores
    """
    
    if Qv is None:
        Q, vv = linalg.eigh( H.dot(H.T) )
    else:
        Q, vv = Qv
        
    # Precomputations. speeds up everything x2.
    Q2 = Q**2
    Qy = Q.T.dot(y)
    
    alphas = [a for a in alphas if (-a) not in vv]


    #@profile
    def Gy(alpha):
        return Q.dot( Qy.T / (vv + alpha))

    #@profile
    def diag_G(alpha):
        v_prime = 1.0 / (vv + alpha)
        return (v_prime * Q2).sum(axis=-1)

    #@profile
    def looe(alpha):
        """ Leave-one-out mean error for alpha """
        errs = Gy(alpha) / diag_G(alpha)
        return (errs**2).mean()

    alphas_scores = []

    if send_state != None:
        percentstep = 10
        nalphastep = len(alphas)/percentstep

    for i,a in enumerate(alphas):
        if (send_state != None) and (i%nalphastep == 0) :
            send_state(min(i/len(alphas)*100,99))
        alphas_scores.append(looe(a))
    
    best_alpha = alphas[np.argmin(alphas_scores)]

    if optimize:
        best_alpha = fmin(looe, best_alpha, disp=0, maxiter=10)[0]

    best_x = H.T.dot( Gy(best_alpha) )

    return best_alpha, best_x, alphas_scores


#@profile
def infer_control(H, y, Nic, alphas=None, eps_L=.0001, positive_solution=False,
                  variances=None, send_state=None):
    """ Infers a control with a smoothness penalty.
    
    
    Returns
    --------
    
    
    u, y_smoothed, ic, best_alpha, alphas_scores
      A set of variables as described below.
    
    
    u
      A vector representing the infered control.
    
    
    y_smoothed
      The predicted value of the observed data at the same time points
      as the data. y_smoothed will appear smoothed compared to y.
    
    
    ic
      initial conditions and other coefficients which are
      not in u.
      
      
    best_alpha
      alpha chosen as a regularization parameters
    
    
    alphas_scores
      Generalized Cross Validation scores for the different
      alphas.
      
    """
        

    if alphas is None:
        alphas = DEFAULT_ALPHAS

    Ny, Nuic = H.shape


    iL = make_iL(Nuic - Nic, Nic, eps=eps_L)
    
    assert iL.sum() != np.nan
    assert H.sum() != np.nan
    
    HiL = H.dot(iL)
    
    if hasattr(alphas, '__iter__'):
        K = HiL.dot(HiL.T)
        vv, Q = linalg.eigh( K )

    
    if positive_solution:
        if not CVXOPT_FOUND:
            raise ValueError ("Install cvxopt for positive solutions. ")
        if hasattr(alphas, '__iter__'):
            alpha, v, scores  = GCV(y, HiL, alphas, Qv = (Q,vv) , send_state=send_state)
        else:
            alpha, scores = alphas, None

        G = cvm( -1.0*np.identity(Nuic))
        h = cvm(1.0*np.zeros((Nuic,1)))
        K = make_L(Nuic - Nic, Nic, eps_L=eps_L)
        result = constrained_ridge( M=H, y=y, alpha=alpha, K=K, G=G, h=h,
                                    invcov=variances)
        ic, uu = result[:Nic], result[Nic:]
        yyhat = H.dot(result)
        return uu, yyhat, ic, alpha, scores
    
    # Non-positive solutions
    #@profile
    def treat_y(y):
        """ GCV on one observation vector """
        if hasattr(alphas, '__iter__'):
            alpha, v, scores  = GCV(y, HiL, alphas, Qv = (Q,vv) , send_state=send_state)
        else:
            alpha, v, scores = alphas, None, None 
        
        K_ = HiL.T.dot(HiL)
        for i in range(len(K_)):
            K_[i,i] += alpha
        iK = linalg.inv(K_)
        iLiK = iL.dot(iK).dot(HiL.T)
        uu = iLiK.dot(y)
        ic, u = uu[:Nic], uu[Nic:]
        yhat = H.dot(uu)
        return u, yhat, ic, alpha, scores
      
    
    if isinstance(y, list):
        return map( treat_y, y)
    else:
        return treat_y(y)


### VARIANTS FOR POSITIVE CURVES


def compute_MM_KK_q(M,y,K=None,invcov=None, K_weigths=None):
    """ Precomputations to minimize |M*x-y|_2 + |K*x|_2 
        subject to ( G*x <= h ) and  ( A*x = b )
        This is a wrapper of cvxopt.
        See documentation of cvxopt.solvers.qp for more. """ 

    Nu= M.shape[1]
    if invcov is None:
        MM = M.T.dot(M)
        q = -M.T.dot(y)
    else:
        if len(invcov.shape) == 1:
            invcov = np.diag(invcov)
        MM = M.T.dot(invcov).dot(M)
        q = -invcov.dot(M).T.dot(y)
        
    if K is None:
        K = np.identity(Nu)
    if K_weigths is None:
        KK = K.T.dot(K)
    else :
        Kw = np.diag(K_weigths)
        KwK = Kw.dot(K)
        KK = KwK.T.dot(KwK)
    
    return MM,KK,q
    
def constrained_ridge(M=None,y=None,alpha=0,K=None,invcov=None,
                      K_weigths=None,G=None, h=None, A=None, b=None,
                      solver=None, initvals=None, precomputed=None):
    """ Minimizes |M*x-y|_2 + |K*x|_2 
        subject to ( G*x <= h ) and  ( A*x = b )
        This is a wrapper of cvxopt.
        See documentation of cvxopt.solvers.qp for more. """ 
    
    if not CVXOPT_FOUND:
      raise ImportError("Constrained inversion requires cvxopt, but I didn't"
                        " find it. Please install cvxopt.")

    if precomputed is not None:
        
        MM,KK,q = precomputed
        
    else:
        
        MM,KK,q = compute_MM_KK_q(M,y,K=K,invcov=invcov,K_weigths=K_weigths)
        
    P = cvm( MM+alpha*KK)
    q.shape = (q.size,1)
    q = cvm(q)
    if G is not None: G = cvm(G)
    if h is not None: h = cvm(h)
    if A is not None: A = cvm(A)
    if b is not None: b = cvm(b)
            
    
    sol = cvxopt.solvers.qp(P,q,G,h,A,b,solver,initvals)
    xx = np.array(sol['x']).flatten()
    return xx


def make_estimator(A,B,C,Ncoefs, tty,yobs,H=None,Lpower = 1,
                   u_func = None, tt_u = None, pwl=True, positive_u=True,
                   invcov=None):
    
    if H is None:
        
        if pwl:
            
            H = make_H_pwl(A,B,C,ttu, tty)
            
        else:
            
            H = make_H(A,B,C,u_func,Ncoefs, tty)
    
    L = make_L(Ncoefs,Lpower)
               
    MM,KK,q = compute_MM_KK_q(M=H,y=yobs,K=L, invcov=invcov) # precomputation
    
    if positive_u:
        
        G = cvm( -np.identity(Ncoefs))
        h = cvm(np.zeros((Ncoefs,1)))
        
        def estimator(alpha):
            
            return constrained_ridge(precomputed=(MM,KK,-q), 
                          alpha=alpha,G=G, h=h, A=None, b=None,
                          invcov=invcov)
    else:
        
        def estimator(alpha):
            
            return alg.solve( MM+alpha*KK,q )
        
    
    return estimator,H,L