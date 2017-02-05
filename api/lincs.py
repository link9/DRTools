def updn_to_str(x):
    if list(set(list(x))) == ['DN', 'UP']:
        return 'UPDN'
    elif list(set(list(x))) == ['UP']:
        return 'UP'
    elif list(set(list(x))) == ['DN']:
        return 'DN'

def json_from_url(url, ret=10):
    
    while True:
        try:
            r = requests.get(url)
            
            if r.status_code == 200:
                return json.loads(r.content)

            else:
                if (ret>1):
                    ret-=1
                    continue

                print('Couldn\'t get json, returning status code')
                return r.status_code
            
        except:
            pass
        



def load_samples(drugname):
    #url = 'http://amp.pharm.mssm.edu/Harmonizome/api/1.0/dataset/LINCS+L1000+CMAP+Signatures+of+Differentially+Expressed+Genes+for+Small+Molecules'
    #content = json_from_url(url)
    
    #use when test mode
    with open('../../data/l1000_samples.json', 'r') as f:
        content = json.loads(f.read())
    
    if type(content) == dict:
    
        geneSets = pd.DataFrame(content['geneSets'])

        geneSets['drug'] = geneSets.name.str.split('_').apply(lambda x : x[1])
        geneSets['sample'] = geneSets.name.str.split('/').apply(lambda x : x[0])
        geneSets = geneSets.loc[geneSets.drug.str.lower() == drugname.lower(),:]
        geneSets.href = geneSets.href.apply(lambda x: 'http://amp.pharm.mssm.edu/Harmonizome'+x)

        geneSets = geneSets.loc[:,['href','drug','sample']]

        return geneSets
    
    else:
        print ('error')
        return -1
    
def sampleToJson(sample):
    url = 'http://amp.pharm.mssm.edu/Harmonizome/api/1.0/gene_set/' + \
            sample + '/LINCS+L1000+CMAP+Signatures+of+Differentially+Expressed+Genes+for+Small+Molecules'
    
    return json_from_url(url)
    
    
    
def drugSamplesToDegs(drugname):
    
    geneSets = load_samples(drugname)  
    samples = []
              
    for idx, val in geneSets['sample'].iteritems():
        samples.append(sampleToDegs(val))
        
    return pd.DataFrame(pd.concat(samples).groupby('gene').updown.apply(updn_to_str))

def geneDrugTable(druglist):
    
    updowns = []
    
    for drugname in druglist:
        updowns.append( drugSamplesToDegs(drugname) )
    
    df_merged = multipleMerge(updowns, how='outer')
    df_merged.columns = druglist
    return df_merged
    
    

def sampleToDegs(sample_str):
    
    print('Fetching : '+ sample_str)
    
    sam = pd.DataFrame(sampleToJson(sample_str)['associations'])  

    sam.thresholdValue = sam.thresholdValue.replace(1.0, "UP")
    sam.thresholdValue = sam.thresholdValue.replace(-1.0, "DN")
    sam.gene = sam.gene.apply(lambda x: x['symbol'])

    sam.columns = ['gene','updown']
    return sam

def multipleMerge(dflist, how):
    df_merged = dflist[0]

    for i in range(1, len(dflist)):
        df_merged = pd.merge(df_merged, dflist[i], how='outer', left_index=True, right_index=True)

    return df_merged
