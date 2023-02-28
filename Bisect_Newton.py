import numpy as np
import matplotlib.pyplot as plt

def bisect_method(f,a,b,tol,nmax,vrb=False):
    count = 0;

    an = a; bn=b; n=0;
    xn = (an+bn)/2;
    
    rn=np.array([xn]);
    r=xn;
    ier=0;

    if vrb:
        print("\n Bisection method with nmax=%d and tol=%1.1e\n" % (nmax, tol));

    
    if f(a)*f(b)>=0:
        print("\n Interval is inadequate, f(a)*f(b)>=0. Try again \n")
        r = "None";
        return r;
    else:
        
        if vrb:
            print("\n|--n--|--an--|--bn--|----xn----|-|bn-an|--|---|f(xn)|---|");

            
            fig, (ax1, ax2) = plt.subplots(1, 2);
            fig.suptitle('Bisection method results');
            ax1.set(xlabel='x',ylabel='y=f(x)');
            
            xl=np.linspace(a,b,100,endpoint=True); yl=f(xl);
            ax1.plot(xl,yl);

        while n<=nmax:
            if vrb:
                print("|--%d--|%1.4f|%1.4f|%1.8f|%1.8f|%1.8f|" % (n,an,bn,xn,bn-an,np.abs(f(xn))));

                
                xint = np.array([an,bn]);
                yint=f(xint);
                ax1.plot(xint,yint,'ko',xn,f(xn),'rs');
                

            
            if (bn-an)<2*tol:
                ier=1;
                break;

            
            if f(an)*f(xn)<0:
                count = count +1;
                bn=xn;
            else:
                count = count +1;
                an=xn;

            
            n += 1;
            xn = (an+bn)/2;
            rn = np.append(rn,xn);

    
    r=xn;

    return [r,count];

def bisect_method_stop(f,fd,a,b,tol,nmax,vrb=False):
    
    count = 0;
    an = a; bn=b; n=0;
    xn = (an+bn)/2;
    
    rn=np.array([xn]);
    r=xn;
    ier=0;

    if vrb:
        print("\n Bisection method with nmax=%d and tol=%1.1e\n" % (nmax, tol));

    
    if f(a)*f(b)>=0:
        print("\n Interval is inadequate, f(a)*f(b)>=0. Try again \n")
        r = "None";
        return r;
    else:
        
        if vrb:
            print("\n|--n--|--an--|--bn--|----xn----|-|bn-an|--|---|f(xn)|---|");

            
            fig, (ax1, ax2) = plt.subplots(1, 2);
            fig.suptitle('Bisection method results');
            ax1.set(xlabel='x',ylabel='y=f(x)');
            
            xl=np.linspace(a,b,100,endpoint=True); yl=f(xl);
            ax1.plot(xl,yl);

        while n<=nmax:
            if vrb:
                print("|--%d--|%1.4f|%1.4f|%1.8f|%1.8f|%1.8f|" % (n,an,bn,xn,bn-an,np.abs(f(xn))));

                
                xint = np.array([an,bn]);
                yint=f(xint);
                ax1.plot(xint,yint,'ko',xn,f(xn),'rs');
                

            
            if (bn-an)<2*tol:
                ier=1;
                break;
            
            if fd(xn) < 1:
                count = count + 1;
                return [xn,count];

            
            if f(an)*f(xn)<0:
                count = count + 1;
                bn=xn;
            else:
                count = count + 1;
                an=xn;

            
            n += 1;
            xn = (an+bn)/2;
            rn = np.append(rn,xn);

    
    r=xn;

    return r;

def newton(f,fp,p0,tol, Nmax):
        count = 0;
        p = np.zeros(Nmax+1);
        p[0] = p0

        for it in range(Nmax):
            count = count + 1;
            p1 = p0-f(p0)/fp(p0)
            p[it+1] = p1
            if (abs(p1-p0) < tol):
                pstar = p1
                info = 0
                return [pstar,count]
            p0 = p1

        pstar = p1
        info = 1

        return [pstar,count]

f = lambda x: np.exp(x**2+7*x-30) - 1
df = lambda x: (2*x + 7)*np.exp(x**2+7*x-30)

[r1,count_b] = bisect_method_stop(f,df,2,4.5,10**-10,200)
[r2,count_n] = newton(f,df,r1,10**-6,200)
count = count_b + count_n

r_b = bisect_method(f,2,4.5,10**-10,200)

r_n = newton(f,df,4.5,10**-10,200)

print("With Bisect Method [root, count]: ", r_b)
print("With Newton Method [root, count]: ", r_n)
print("With combined method [root, count]: ", [r2,count])