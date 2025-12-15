import sauna

# From COMMARA
covariances = sauna.Covariances()
covariances.library = 'COMMARA'
covariances.group_structure = sauna.GROUP_STRUCTURES['ECCO-33']
covariances.from_commara('../../covariances/COMMARA.cov')

print(covariances)
print(covariances.covariances[0])