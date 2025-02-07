
import torch
from cte import L


def Fz(z):
    z_tensor = torch.as_tensor(z)
    L_tensor = torch.as_tensor(L)

    condition = (z_tensor >= L_tensor / 12) & (z_tensor <= 11 * L_tensor / 12)
    expoente = -5 * ((z_tensor - (0.5 * L_tensor)) ** 2 * 1 / L_tensor ** 2)

    value = torch.where(condition,
                        torch.exp(expoente) * torch.cos(6 * torch.pi * z_tensor / L_tensor),
                        torch.tensor(0))
    return value

def F_z1(z, L):
    z_tensor = torch.as_tensor(z)
    L_tensor = torch.as_tensor(L)
    
    condition = (z_tensor >= 3*L_tensor/8) & (z_tensor <= 5*L_tensor/8)

    value = torch.where(condition, 
                        torch.tensor(1),
                        torch.tensor(0))
    return(value)

def F_z2(z):

    z_tensor = torch.as_tensor(z)
    L_tensor = torch.as_tensor(L)

    condition = (z_tensor >= 2*L_tensor/8) & ( z_tensor <= 6*L_tensor/8)
    expoente = (10**(4))*(z_tensor - 6*L_tensor/8)/2
    
    value = torch.where(condition,
                        torch.exp(expoente),
                        torch.tensor(0))
    
    return(value)


def F_z3(z):

    z_tensor = torch.as_tensor(z)
    L_tensor = torch.as_tensor(L)

    condition = (z_tensor >= L_tensor / 12) & (z_tensor <= 11 * L_tensor / 12)
    expoente = -5 * ((z_tensor - (0.5 * L_tensor)) ** 2 * 1 / L_tensor ** 2)

    value = torch.where(condition,
                        torch.exp(expoente) * torch.cos(6 * torch.pi * z_tensor / L_tensor),
                        torch.tensor(0))
    return value


def F_z4(z):
    z_tensor = torch.as_tensor(z)
    L_tensor = torch.as_tensor(L)

    condition = ((z_tensor >= 1/4*L_tensor) & ( z_tensor <= 3/4*L_tensor))
    
    value = torch.where(condition,
                        torch.tensor((torch.special.bessel_j0((1.6)*10**(-6) * z_tensor))),
                        torch.tensor(0))
    
    return(value)


#funcao artigo 
def F_z5(z):

    z_tensor = torch.as_tensor(z)
    L_tensor = torch.as_tensor(L)

    l1 = (1/5*L_tensor - L_tensor/50)
    l2 = (1/5*L_tensor + L_tensor/50)

    l3 = (1/2*L_tensor - L_tensor/10)
    l4 = (1/2*L_tensor + L_tensor/10)

    l5 = (4/5*L_tensor - L_tensor/50)
    l6 = (4/5*L_tensor + L_tensor/50)


    condition = ((z_tensor >= l1) & ( z_tensor <= l2))
    
    condition1 = ((z_tensor >= l1) & (z_tensor <= l2))
    condition2 = ((z_tensor >= l3) & (z_tensor <= l4))
    condition3 = ((z_tensor >= l5) & (z_tensor <= l6))


    value = torch.where(condition1,
                        torch.tensor((-4 * (z_tensor - l1) * (z_tensor - l2)) / (l2 - l1)**2),
                        torch.where(condition2, 1,
                                    torch.where(condition3,
                                                torch.tensor((-4 * (z_tensor - l5) * (z_tensor - l6)) / (l6 - l5)**2), 0)
                                    )
                        )

        
    return(value)
    
def F_z6(z):

    z_tensor = torch.as_tensor(z)
    
    condition1 = (z_tensor >= torch.tensor(10*10**(-2))) & (z_tensor <= torch.tensor(15*10**(-2)))
    condition2 = (z_tensor >= torch.tensor(20*10**(-2))) & (z_tensor <= torch.tensor(25*10**(-2)))
    

    value = torch.where(condition1 , torch.tensor(1), 
                        torch.where(condition2, torch.tensor(1), torch.tensor(0)))

    return(value)