import sauna

covariances = sauna.Covariances()
covariances.library = 'ENDF-B-VII.1'
covariances.group_structure = sauna.GROUP_STRUCTURES['ABBN-28']
covariances.from_excels( f'../../covariances/{covariances.library}-28/')

print(covariances)
print(covariances.covariances[0])