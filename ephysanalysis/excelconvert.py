import dill as pickle
import pandas
pfile = open('VGATDCNMaps.p','rb')
data=pickle.load(pfile)
pfile.close()
writefile=pandas.ExcelWriter('VGATDCNMaps.xlsx')
data.to_excel(writefile)
writefile.save()
