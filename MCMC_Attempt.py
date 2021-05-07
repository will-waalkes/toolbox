# -*- coding: utf-8 -*-
"""
Created on Thu May 18 13:23:23 2017

@author: willwaalkes
"""
import os
import sys
import numpy as np
import scipy as sp
import pandas as pd #This package lets me read in files with row or column labels!
import matplotlib.pyplot as plt
import seaborn as sns #Based on matplot lib, but allows different statistical graphs to be visualized

from scipy.stats import norm

sns.set_style('white')
sns.set_context('talk')

np.random.seed(123)

'''
Lets generate some data: 100 points from a normal centered around zero.
Our goal will be to estimate the posterior of the mean mu
(we'll assume that we know the standard deviation to be 1).
'''

data = np.random.randn(20)

#ax = plt.subplot()
#sns.distplot(data, kde=True, ax=ax) #kde adds a best fit line!
#_ = ax.set(title='Histogram of observed data', xlabel='x', ylabel='# observations');


'''
We have created a data set. Next, we have to define our model.
In this simple case, we will assume that this data is normal distributed,
i.e. the likelihood of the model is normal. As you know, a normal distribution has two parameters
-- mean $\mu$ and standard deviation $\sigma$.
For simplicity, we'll assume we know that $\sigma = 1$
and we'll want to infer the posterior for $\mu$. For each parameter we want to infer,
we have to chose a prior. For simplicity, lets also assume a Normal distribution as a prior for $\mu$.

What is convenient, is that for this model, we actually can compute the posterior analytically.
That's because for a normal likelihood with known standard deviation, the normal prior for mu is
conjugate (conjugate here means that our posterior will follow the same distribution as the prior),
so we know that our posterior for $\mu$ is also normal. We can easily look up on wikipedia
how we can compute the parameters of the posterior.


'''
#print data
#print data.sum()
#os.sys.exit()

def calc_posterior_analytical(data, x, mu_0, sigma_0):
    sigma = 1.
    n = len(data)
    mu_post = (mu_0 / sigma_0**2 + data.sum() / sigma**2) / (1. / sigma_0**2 + n / sigma**2)
    sigma_post = (1. / sigma_0**2 + n / sigma**2)**-1
    return norm(mu_post, np.sqrt(sigma_post)).pdf(x)

#ax = plt.subplot()
#x = np.linspace(-1, 1, 500) #The visible window of the probability distribution function
#posterior_analytical = calc_posterior_analytical(data, x, 0., 1.)
#ax.plot(x, posterior_analytical)
#ax.set(xlabel='mu', ylabel='belief', title='Analytical posterior');
#sns.despine() #all borders or just some

'''
This shows our quantity of interest, the probability of $\mu$'s values after having
seen the data, taking our prior information into account. Lets assume, however, that
our prior wasn't conjugate and we couldn't solve this by hand which is usually the case.

Explaining MCMC sampling with code

Now on to the sampling logic. At first, you find starting parameter position
(can be randomly chosen), lets fix it arbitrarily to:
'''
#mu_current = 1.
'''
Then, you propose to move (jump) from that position somewhere else (that's the Markov part).
You can be very dumb or very sophisticated about how you come up with that proposal.
The Metropolis sampler is very dumb and just takes a sample from a normal distribution
(no relationship to the normal we assume for the model) centered around your
current mu value (i.e. mu_current) with a certain standard deviation (proposal_width)
that will determine how far you propose jumps (here we're use scipy.stats.norm):
'''
#proposal = norm(mu_current, proposal_width).rvs()

'''
Next, you evaluate whether that's a good place to jump to or not. If the resulting
normal distribution with that proposed mu explaines the data better than your old mu,
you'll definitely want to go there. What does "explains the data better" mean?
We quantify fit by computing the probability of the data, given the likelihood (normal)
with the proposed parameter values (proposed mu and a fixed sigma = 1). This can easily
be computed by calculating the probability for each data point using
scipy.stats.normal(mu, sigma).pdf(data) and then multiplying the individual probabilities,
i.e. compute the likelihood (usually you would use log probabilities but we omit this here):
'''

#likelihood_current = norm(mu_current, 1).pdf(data).prod()
#likelihood_proposal = norm(mu_proposal, 1).pdf(data).prod()

# Compute prior probability of current and proposed mu        
#prior_current = norm(mu_prior_mu, mu_prior_sd).pdf(mu_current)
#prior_proposal = norm(mu_prior_mu, mu_prior_sd).pdf(mu_proposal)

# Numerator of Bayes formula
#p_current = likelihood_current * prior_current
#p_proposal = likelihood_proposal * prior_proposal

'''
Up until now, we essentially have a hill-climbing algorithm that would just propose movements
into random directions and only accept a jump if the mu_proposal has higher likelihood than
mu_current. Eventually we'll get to mu = 0 (or close to it) from where no more moves will
be possible. However, we want to get a posterior so we'll also have to sometimes accept
moves into the other direction. The key trick is by dividing the two probabilities,
'''
#p_accept = p_proposal / p_current
'''we get an acceptance probability. You can already see that if p_proposal is larger,
that probability will be > 1 and we'll definitely accept. However, if p_current is larger,
say twice as large, there'll be a 50% chance of moving there:
'''
#accept = np.random.rand() < p_accept

#if accept:
    # Update position
#    cur_pos = proposal
#This simple procedure gives us samples from the posterior, and P(x) gets cancelled out!

def sampler(data, samples=4, mu_init=.5, proposal_width=.5, plot=True, mu_prior_mu=0, mu_prior_sd=1.):
    mu_current = mu_init
    posterior = [mu_current]
    for i in range(samples):
        # suggest new position
        mu_proposal = norm(mu_current, proposal_width).rvs()

        # Compute likelihood by multiplying probabilities of each data point
        likelihood_current = norm(mu_current, 1).pdf(data).prod()
        likelihood_proposal = norm(mu_proposal, 1).pdf(data).prod()
        
        # Compute prior probability of current and proposed mu        
        prior_current = norm(mu_prior_mu, mu_prior_sd).pdf(mu_current)
        prior_proposal = norm(mu_prior_mu, mu_prior_sd).pdf(mu_proposal)
        
        p_current = likelihood_current * prior_current
        p_proposal = likelihood_proposal * prior_proposal
        
        # Accept proposal?
        p_accept = p_proposal / p_current
        
        # Usually would include prior probability, which we neglect here for simplicity
        accept = np.random.rand() < p_accept
        
        if plot:
            plot_proposal(mu_current, mu_proposal, mu_prior_mu, mu_prior_sd, data, accept, posterior, i)
        
        if accept:
            # Update position
            mu_current = mu_proposal
        
        posterior.append(mu_current)
        
    return posterior

# Function to display
def plot_proposal(mu_current, mu_proposal, mu_prior_mu, mu_prior_sd, data, accepted, trace, i):
    from copy import copy
    trace = copy(trace)
    fig, (ax1, ax2, ax3, ax4) = plt.subplots(ncols=4, figsize=(16, 4))
    fig.suptitle('Iteration %i' % (i + 1))
    x = np.linspace(-3, 3, 5000)
    color = 'g' if accepted else 'r'
        
    # Plot prior
    prior_current = norm(mu_prior_mu, mu_prior_sd).pdf(mu_current)
    prior_proposal = norm(mu_prior_mu, mu_prior_sd).pdf(mu_proposal)
    prior = norm(mu_prior_mu, mu_prior_sd).pdf(x)
    ax1.plot(x, prior)
    ax1.plot([mu_current] * 2, [0, prior_current], marker='o', color='b')
    ax1.plot([mu_proposal] * 2, [0, prior_proposal], marker='o', color=color)
    ax1.annotate("", xy=(mu_proposal, 0.2), xytext=(mu_current, 0.2),
                 arrowprops=dict(arrowstyle="->", lw=2.))
    ax1.set(ylabel='Probability Density', title='current: prior(mu=%.2f) = %.2f\nproposal: prior(mu=%.2f) = %.2f' % (mu_current, prior_current, mu_proposal, prior_proposal))
    
    # Likelihood
    likelihood_current = norm(mu_current, 1).pdf(data).prod()
    likelihood_proposal = norm(mu_proposal, 1).pdf(data).prod()
    y = norm(loc=mu_proposal, scale=1).pdf(x)
    sns.distplot(data, kde=False, norm_hist=True, ax=ax2)
    ax2.plot(x, y, color=color)
    ax2.axvline(mu_current, color='b', linestyle='--', label='mu_current')
    ax2.axvline(mu_proposal, color=color, linestyle='--', label='mu_proposal')
    #ax2.title('Proposal {}'.format('accepted' if accepted else 'rejected'))
    ax2.annotate("", xy=(mu_proposal, 0.2), xytext=(mu_current, 0.2),
                 arrowprops=dict(arrowstyle="->", lw=2.))
    ax2.set(title='likelihood(mu=%.2f) = %.2f\nlikelihood(mu=%.2f) = %.2f' % (mu_current, 1e14*likelihood_current, mu_proposal, 1e14*likelihood_proposal))
    
    # Posterior
    posterior_analytical = calc_posterior_analytical(data, x, mu_prior_mu, mu_prior_sd)
    ax3.plot(x, posterior_analytical)
    posterior_current = calc_posterior_analytical(data, mu_current, mu_prior_mu, mu_prior_sd)
    posterior_proposal = calc_posterior_analytical(data, mu_proposal, mu_prior_mu, mu_prior_sd)
    ax3.plot([mu_current] * 2, [0, posterior_current], marker='o', color='b')
    ax3.plot([mu_proposal] * 2, [0, posterior_proposal], marker='o', color=color)
    ax3.annotate("", xy=(mu_proposal, 0.2), xytext=(mu_current, 0.2),
                 arrowprops=dict(arrowstyle="->", lw=2.))
    #x3.set(title=r'prior x likelihood $\propto$ posterior')
    ax3.set(title='posterior(mu=%.2f) = %.5f\nposterior(mu=%.2f) = %.5f' % (mu_current, posterior_current, mu_proposal, posterior_proposal))
    
    if accepted:
        trace.append(mu_proposal)
    else:
        trace.append(mu_current)
    ax4.plot(trace)
    ax4.set(xlabel='iteration', ylabel='mu', title='trace')
    plt.tight_layout()
    #plt.legend()



