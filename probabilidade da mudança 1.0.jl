

sigma_g = 0.02
sigma_c = 0.10
mu= 0.10
sigma=0.05
sigma_1=0.10
T_=20.0
tau=10.0
t=tau -1.0
gamma = 5.0
sigma_gl = 0.01
sigma_gh = 0.03
mu_gl =-0.008
mu_gh =0.008
cl = -(sigma_c^2)/2
ch=cl
c_til = cl/((gamma-1.0)*(T_-tau))
mu_til_gh = mu_gh - sigma_gh^2*(T_-tau)*(gamma - 1)/2
mu_til_gl = mu_gl - sigma_gl^2*(T_-tau)*(gamma - 1)/2
sigma_til= sigma_g/((gamma-1)*(T_-tau))
sigma_til_gh = sigma_gh/((gamma-1)*(T_-tau))
sigma_t = 1/((1/sigma_g^2) +1/(sigma^2))
#sigma_ct= 1/((1/sig) )


using Distributions
g_hat=LinRange(-0.03,0.03,300)
ph=zeros(300)
for i in eachindex(g_hat)
   # g_hat=-0.03
   a = cdf(Normal(c_til , sigma_til_gh),c_til+mu_til_gl- mu_til_gh)
b = cdf(Normal(g_hat[i] - sigma_gh^2*((gamma-1)*(t-tau))/2 , sqrt(sigma_g^2 -sigma_t^2)), mu_til_gh - c_til)

c = 1/(sqrt(2*pi)*sigma_til_gh)

using QuadGK
k=integral, err = quadgk(x -> (1-a)*b*c *exp(-( 0.5*((x-c_til)/sigma_til_gh)^2)), -Inf, +Inf, rtol=1e-8)
k
ph[i]=k[1]

a = cdf(Normal(c_til , sigma_til_gl),c_til+mu_til_gh- mu_til_gl)
b = cdf(Normal(g_hat[i] - sigma_gl^2*((gamma-1)*(t-tau))/2 , sqrt(sigma_g^2 -sigma_t^2)), mu_til_gl - c_til)

c = 1/(sqrt(2*pi)*sigma_til_gl)
k=integral, err = quadgk(x -> (1-a)*b*c *exp(-( 0.5*((x-c_til)/sigma_til_gl)^2)), -Inf, +Inf, rtol=1e-8)
pl[i]=k[1]
end
ph

using Plots
plot(g_hat, ph)
?plot


