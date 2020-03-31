syms m a g theta(t)
eqn = m*a == -m*g*sin(theta)

syms r
eqn = subs(eqn,a,r*diff(theta,2))

eqn = isolate(eqn,diff(theta,2))

eqn