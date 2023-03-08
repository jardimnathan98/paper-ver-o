
sigma_g =0.02
sigma_c = 0.10
mu= 0.10
sigma=0.05
sigma_1=0.10
t=20
tau=9
gamma = 5
sigma_gl = 0.01
sigma_gh = 0.03
mu_gl =-0.008
mu_gh =0.008
cl = -(sigma_c^2)/2 
ch=cl
c_til = cl/(gamma-1)*(t-tau)
mu_til_gh = mu_gh - sigma_gh^2*(t-tau)*(gamma - 1)/2
mu_til_gl = mu_gl - sigma_gl^2*(t-tau)*(gamma - 1)/2
sigma_til= sigma_g/((gamma-1)*(t-tau))
sigma_til_gh= sigma_gh/((gamma-1)*(t-tau))

 a<- pnorm(c_til+mu_til_gl- mu_til_gh , mean = c_til, sd = sigma_til_gh) 
g_hat=0.0
b<-  pnorm( mu_til_gh - c_til , mean = (g_hat - sigma_gh^2*((gamma-1)*(t-tau))/2), sd = sqrt(sigma_1^2 - sigma_gh^2))


c= 1/(sqrt(2*pi)*sigma_til_gh)

integrand <- function(x) {(1-a)*b*c *exp(-( 0.5*((x-c_til)/sigma_til_gh)^2)) }
integrate(integrand, lower = -Inf, upper = Inf)
