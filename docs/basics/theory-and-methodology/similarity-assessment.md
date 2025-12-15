---
icon: reflect-horizontal
---

# Similarity Assessment

Choosing experiments for validation is widely based upon the similarity indices such as $$c_k$$, $$E$$, and $$G$$. They allow numerically quantifying how much two systems are similar and finding the amount of shared information. The most accepted index is $$c_k$$, the correlation coefficient, defined as:

$$
c_k=\frac{S_aCS^T_e}{\sqrt{S_aCS_a^T}\sqrt{S_eCS_e^T}}
$$

where $$S_a$$ is the application sensitivity vector;  $$S_e$$ is the application sensitivity vector; $$C$$ is the covariance matrix.

The second index is $$E$$ and defined as the cosine between two vectors:

$$
E=\frac{S_aS^T_e}{|S_a||S_e|}
$$

The $$G$$ similarity index is intended to quantify the coverage of the application by an experiment and defined as follows:

$$
G=1-\frac{\sum_{n,x,g}|S_{a,n,x,g}-S_{e',n,x,g}|}{\sum_{n,x,g}|S_{a,n,x,g}|}
$$

where $$S_{a,n,x,g}$$ is the application sensitivity to group $$g$$ of reaction $$x$$ of nuclide $$n$$; $$S_{a,n,x,g}$$ is the experimental sensitivity to group $$g$$ of reaction $$x$$ of nuclide $$n$$;

&#x20;$$S_{e',n,x,g}=S_{e,n,x,g}$$ if $$|S_{a,n,x,g}|\geq |S_{e,n,x,g}|$$ and $$S_{a,n,x,g}/|S_{a,n,x,g}|=S_{e,n,x,g}/|S_{e,n,x,g}|$$;&#x20;

&#x20;$$S_{e',n,x,g}=S_{a,n,x,g}$$ if $$|S_{a,n,x,g}|\lt |S_{e,n,x,g}|$$ and $$S_{a,n,x,g}/|S_{a,n,x,g}|=S_{e,n,x,g}/|S_{e,n,x,g}|$$;

$$S_{e',n,x,g}=0$$ otherwise.

