from flask import Flask
import psycopg2
import psycopg2.extras
import argparse
import json
import numpy as np
import datetime
import zlib

conn,cur = None,None

resolutions = ['microrregiao','mesorregiao','municipio','estado']

classes = ['verylow','low','mid','high','veryhigh']

clippings = ['semiarido']

colorMap = ['#d7191c', '#fdae61', '#ffffbf', '#a6d96a', '#1a9641']

mapTemplates = {}

cache = {}
cacheType = -1

currentYear = datetime.datetime.today().year

def startApp():
    app = Flask('sismoiWS')
    return app

app = startApp()

def ProcessCmdLine():
    parser = argparse.ArgumentParser(description="SISMOI WebServices.")
    parser.add_argument("-d", "--debug", help="Activate Flask debug mode", action='store_true')
    parser.add_argument("-host", "--host", help="Host IP", type=str, default="127.0.0.1")
    parser.add_argument("-p", "--port", type=str, help="Port to be used", default=5000)
    parser.add_argument("-c", "--cachetype", type=int, help="Cache type: 0 for no cache, 1 for plain, 2 for compressed",
                        default=1)
    return parser.parse_args()

def connect():
    conn = psycopg2.connect("dbname='sismoi' user='sismoi' host='200.133.39.41' password='142857'")
    cur = conn.cursor()
    return conn, cur

def getValue(sql):
    cur.execute(sql)
    row=cur.fetchone()
    if row != None:
        return row[0]
    else:
        return None

def findElement(data,keyfield, feature):
    for rec in data:
        if rec[keyfield] == int(feature['properties'][keyfield]):
            return rec
    raise Exception('Element not found in findElement. keyfield: {0} Feature: \n {1}'.format(keyfield,feature))

def getDictResultset(sql):
    dictcur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    dictcur.execute(sql)
    resultset =dictcur.fetchall()
    dictResult = []
    for row in resultset:
        dictResult.append(dict(row))
    return dictResult

def getStates():
    cur.execute("select distinct state from county where state is not null")
    return([token[0] for token in cur.fetchall()])

def getYears():
    cur.execute("""select indicator_id,string_agg(year::character varying,',') as years from 
                   (select distinct indicator_id,year from value) a
                   group by indicator_id
                   order by indicator_id
                """)
    return {id: years for id, years in cur.fetchall()}

def toCache(label,value):
    if cacheType == 1:
        cache[label] = value
    elif cacheType == 2:
        cache[label] = zlib.compress(bytes(value,'utf-8'))

def fromCache(label):
    if cacheType == 1:
        return cache[label]
    elif cacheType == 2:
        return zlib.decompress(cache[label])

def inCache(label):
    return label in cache

def featureColor(value, pessimist):
    value = value if (not pessimist) else 1 - value
    return colorMap[0] if (value >= 0.0 and value <= 0.2) else\
           colorMap[1] if (value >  0.2 and value <= 0.4) else\
           colorMap[2] if (value >  0.4 and value <= 0.6) else\
           colorMap[3] if (value >  0.6 and value <= 0.8) else\
           colorMap[4]

def validateClippingResolution(sparams):
    errorMsg=''
    params = dict(token.split('=') for token in sparams.split(','))
    if not(params['resolution'] in resolutions):
        errorMsg='SISMOI Err: Resolução errada: {0}.'.format(params['resolution'])
    if params['clipping'] not in clippings:
        errorMsg="É necessário especificar clipping como 'semiarido' ou uma UF válida."
    return errorMsg,params

def validateParams(sparams):
    errorMsg,params=validateClippingResolution(sparams)
    if errorMsg != '':
        return errorMsg, params
    if not 'indicator_id' in params:
        errorMsg = "SISMOI Err: É obrigatório especificar o indicador (indicator_id)."
    elif not 'scenario_id' in params:
        errorMsg = "SISMOI Err: É obrigatório especificar o cenário (scenario_id)."
    elif not params['scenario_id'] in ['1','2','null']:
        errorMsg='Valor de scenario_id inválido, deve ser 1, 2 ou null.'
    elif ('year' in params):
        if not params['year'] in years[int(params['indicator_id'])].split(','):
            errorMsg='SISMOI Err: Ano {0} não existe para o indicador {1}. Anos válidos: {2}'.format(params['year'],
                                                                                       params['indicator_id'],
                                                                                       years[int(params['indicator_id'])])
        elif (params['scenario_id'] == 'null') and (int(params['year']) > currentYear):
            errorMsg='SISMOI Err: É necessário especificar o cenário para ano maior que {0}.'.format(currentYear)
        elif (params['scenario_id'] != 'null') and (int(params['year']) <= currentYear):
            errorMsg='SISMOI Err: Cenários não podem ser especificados para anos menores ou iguais a {0}.'.format(currentYear)
    elif params['scenario_id'] == 'null':
        errorMsg='SISMOI Err: Se o ano não foi especificado, o cenário não pode ser nulo.'
    return errorMsg,params

def addFeatureColor(data):
    for rec in data:
        rec['valuecolor'] = featureColor(rec['value'],rec['pessimist'])
    return data

def toGroupedDict(data,pessimist):
    ret = {}
    for i in range(0,len(data)):
        if not(data[i]['year'] in ret.keys()):
            ret[data[i]['year']] = {}
            for _class,colorFake in zip(classes,np.arange(0.1,1.0,0.2)):
                ret[data[i]['year']].update({_class: {'data': [], 'count': 0,
                                                      'valuecolor': featureColor(colorFake,pessimist)},'count': 0})
        if 'id' in data[i]:
            ret[data[i]['year']][data[i]['class']]['data'].append({'id':data[i]['id'],
                'name':data[i]['name'], 'value':data[i]['value']})
        else:
            ret[data[i]['year']][data[i]['class']]['data'].append({'state':data[i]['state'],
                'value':data[i]['value']})
        ret[data[i]['year']][data[i]['class']]['count']+=1
        ret[data[i]['year']]['count']+=1
    return ret

def getIndicatorByCounty(params):
    data=getDictResultset('''select v.county_id as id, v.indicator_id, v.year, v.scenario_id, i.pessimist, v.value from value v 
                             inner join county c
                             on v.county_id = c.id
                             inner join indicator i
                             on v.indicator_id = i.id
                             where v.indicator_id = {0}
                               and v.scenario_id {1} 
                               {2}
                               and year = {3}
                             order by v.county_id, v.year, v.scenario_id
                          '''.format(params['indicator_id'],
                                     ' = ' +params['scenario_id'] if params['scenario_id'] != 'null' else 'is null',
                                     "and c.state = '{0}'".format(params['clipping']) if params['clipping'] != 'semiarido' else '',
                                     params['year'])
                          )
    return addFeatureColor(data)

def getIndicatorByMicroregion(params):
    data=getDictResultset('''select c.microregion_id as id, m.name, indicator_id, i.pessimist, v.scenario_id,year, avg(v.value) as value
                              from value v
                             inner join county c
                                on v.county_id = c.id 
                             inner join indicator i
                             on v.indicator_id = i.id
                             inner join microregion m
                                on c.microregion_id = m.id
                             where indicator_id = {0}
                               and scenario_id {1}
                               {2}
                               and year = {3}
                            group by  c.microregion_id, m.name, v.indicator_id, i.pessimist, v.year, v.scenario_id
                             order by                   m.name, v.year, v.scenario_id
                          '''.format(params['indicator_id'],
                                     ' = ' +params['scenario_id'] if params['scenario_id'] != 'null' else 'is null',
                                     "and c.state = '{0}'".format(params['clipping']) if params['clipping'] != 'semiarido' else '',
                                     params['year'])
                          )
    return addFeatureColor(data)

def getIndicatorByMesoregion(params):
    data=getDictResultset('''select c.mesoregion_id as id, m.name, indicator_id, pessimist, v.scenario_id, year, avg(v.value) as value
                               from value v
                              inner join indicator i
                                 on v.indicator_id = i.id
                              inner join county c
                                 on v.county_id = c.id 
                              inner join mesoregion m
                                 on c.mesoregion_id = m.id
                              where indicator_id = {0}
                                and scenario_id {1}
                                {2}
                                and year = {3}
                            group by  c.mesoregion_id, m.name, v.indicator_id, i.pessimist, v.year, v.scenario_id
                            order by                   m.name, v.year, v.scenario_id
                          '''.format(params['indicator_id'],
                                     ' = ' +params['scenario_id'] if params['scenario_id'] != 'null' else 'is null',
                                     "and c.state = '{0}'".format(params['clipping']) if params['clipping'] != 'semiarido' else '',
                                     params['year'])
                          )
    return addFeatureColor(data)

def getIndicatorByState(params):
    data=getDictResultset('''select c.state,v.indicator_id, i.pessimist, v.scenario_id, year, avg(v.value) as value
                               from value v
                              inner join indicator i
                                 on v.indicator_id = i.id
                              inner join county c
                                 on v.county_id = c.id 
                              inner join mesoregion m
                                 on c.mesoregion_id = m.id
                              where indicator_id = {0}
                                and scenario_id {1}
                                {2}
                                and year = {3}
                            group by  c.state, v.indicator_id, i.pessimist, v.scenario_id,year
                             order by c.state, v.indicator_id, v.scenario_id,year
                          '''.format(params['indicator_id'],
                                     ' = ' +params['scenario_id'] if params['scenario_id'] != 'null' else 'is null',
                                     "and c.state = '{0}'".format(params['clipping']) if params['clipping'] != 'semiarido' else '',
                                     params['year'])
                          )
    return addFeatureColor(data)

def getTotalByState(params):
    pessimist = getValue('select pessimist from indicator where id = {0}'.format(params['indicator_id']))
    rawdata=getDictResultset('''select c.state, year, 
                             case
                                when avg(v.value) between 0   and 0.2 then 'verylow'
                                when avg(v.value) between 0.2 and 0.4 then 'low'
                                when avg(v.value) between 0.4 and 0.6 then 'mid'
                                when avg(v.value) between 0.6 and 0.8 then 'high'
                                when avg(v.value) > 0.8               then 'veryhigh'
                             end as class, avg(value) as value
                                from 
                                value v
                                inner join county c
                                on v.county_id = c.id
                                inner join indicator i
                                on v.indicator_id = i.id
                                inner join mesoregion m
                                on c.mesoregion_id = m.id
                                where indicator_id = {0}
                                {1}
                                {2} 
                                {3}                
                             group by c.state, year, i.pessimist
                                order by year {4}, c.state
                          '''.format(params['indicator_id'],
                             (' and (scenario_id = {0} or scenario_id is null)'.format(params['scenario_id'])) if params[ 'scenario_id'] != 'null' else '',
                              "and c.state = '{0}'".format(params['clipping']) if params['clipping'] != 'semiarido' else '',
                              'and v.year = {0}'.format(params['year']) if 'year' in params else '',
                              ('desc' if pessimist == 0 else '')
                             )
                          )
    data=toGroupedDict(rawdata,pessimist)
    return data

def getTotalByMesoregion(params):
    pessimist = getValue('select pessimist from indicator where id = {0}'.format(params['indicator_id']))
    rawdata=getDictResultset('''select c.mesoregion_id as id, m.name, year,
                             case
                                when avg(v.value) between 0   and 0.2 then 'verylow'
                                when avg(v.value) between 0.2 and 0.4 then 'low'
                                when avg(v.value) between 0.4 and 0.6 then 'mid'
                                when avg(v.value) between 0.6 and 0.8 then 'high'
                                when avg(v.value) > 0.8               then 'veryhigh'
                             end as class, avg(value) as value
                                from value v
                                inner join county c
                                on v.county_id = c.id
                                inner join indicator i
                                on v.indicator_id = i.id
                                inner join mesoregion m
                                on c.mesoregion_id = m.id
                                where indicator_id = {0}
                                {1}
                                {2}
                                {3}
                             group by year, mesoregion_id, m.name
                                order by year, value {4}, name
                          '''.format(params['indicator_id'],
                             (' and (scenario_id = {0} or scenario_id is null)'.format(params['scenario_id'])) if params[ 'scenario_id'] != 'null' else '',
                              "and c.state = '{0}'".format(params['clipping']) if params['clipping'] != 'semiarido' else '',
                              'and v.year = {0}'.format(params['year']) if 'year' in params else '',
                              ('desc' if pessimist == 0 else '')
                             )
                          )
    data=toGroupedDict(rawdata,pessimist)
    return data

def getTotalByMicroregion(params):
    pessimist = getValue('select pessimist from indicator where id = {0}'.format(params['indicator_id']))
    rawdata=getDictResultset('''select c.microregion_id as id, m.name, year, 
                             case
                                when avg(v.value) between 0   and 0.2 then 'verylow'
                                when avg(v.value) between 0.2 and 0.4 then 'low'
                                when avg(v.value) between 0.4 and 0.6 then 'mid'
                                when avg(v.value) between 0.6 and 0.8 then 'high'
                                when avg(v.value) > 0.8               then 'veryhigh'
                             end as class, avg(value) as value
                                from value v
                                inner join county c
                                on v.county_id = c.id
                                inner join indicator i
                                on v.indicator_id = i.id
                                inner join microregion m
                                on c.microregion_id = m.id
                                where indicator_id = {0}
                                {1}
                                {2}
                                {3}
                             group by year, c.microregion_id, m.name
                                order by year, value {4}, name
                          '''.format(params['indicator_id'],
                             (' and (scenario_id = {0} or scenario_id is null)'.format(params['scenario_id'])) if params[ 'scenario_id'] != 'null' else '',
                              "and c.state = '{0}'".format(params['clipping']) if params['clipping'] != 'semiarido' else '',
                              'and v.year = {0}'.format(params['year']) if 'year' in params else '',
                              ('desc' if pessimist == 0 else '')
                             )
                          )
    data=toGroupedDict(rawdata,pessimist)
    return data

def getTotalByCounty(params):
    pessimist = getValue('select pessimist from indicator where id = {0}'.format(params['indicator_id']))
    rawdata=getDictResultset('''select c.id, c.name, year, v.value,
                             case
                                when value between 0   and 0.2 then 'verylow'
                                when value between 0.2 and 0.4 then 'low'
                                when value between 0.4 and 0.6 then 'mid'
                                when value between 0.6 and 0.8 then 'high'
                                when value > 0.8               then 'veryhigh'
                             end as class
                                from value v
                                inner join county c
                                on v.county_id = c.id
                                inner join indicator i
                                on v.indicator_id = i.id
                                where indicator_id = {0}
                                {1}
                                {2}
                                {3}
                                order by year, value {4}, name
                          '''.format(params['indicator_id'],
                             (' and (scenario_id = {0} or scenario_id is null)'.format(params['scenario_id'])) if params[ 'scenario_id'] != 'null' else '',
                              "and c.state = '{0}'".format(params['clipping']) if params['clipping'] != 'semiarido' else '',
                              'and v.year = {0}'.format(params['year']) if 'year' in params else '',
                              ('desc' if pessimist == 0 else '')
                             )
                          )
    data=toGroupedDict(rawdata,pessimist)
    return data

@app.route("/sismoi/getCacheSize", methods=['GET'])
def getCacheSize():
    return "Not implemented"

@app.route("/sismoi/clearCache", methods=['GET'])
def clearCache():
    global cache
    cache = {}
    return 'New cache size: {0}'.format(getCacheSize())

@app.route("/sismoi/getHierarchy", methods=['GET'])
def getHierarchy():
    if inCache('getHierarchy@'):
        return fromCache('getHierarchy@')
    try:
        data=getDictResultset("""select a.id,a.name,a.title,a.shortname,
                                        a.simple_description,a.complete_description,
                                        a.equation,a.level,a.pessimist,
                                        string_agg(distinct b.indicator_id_master::character varying,',') as indicator_id_master, 
                                        string_agg(year::character varying,',') as years 
                                        from indicator a
                                 left join indicator_indicator b
                                        on b.indicator_id_detail = a.id
                                 left join (select distinct indicator_id,year from value) v
                                        on a.id = v.indicator_id
                                     where level = 1 or indicator_id_master is not null
                                  group by a.id,a.name,a.title,a.shortname,a.simple_description,a.complete_description,
                                           a.equation,a.level,a.pessimist
                                  order by level,id
                              """)
        ret = json.dumps(data)
        toCache('getHierarchy@',ret)
        return ret
    except Exception as e:
        return str(e)

@app.route('/sismoi/getGeometry/<sparams>', methods=['GET'])
def getGeometry(sparams):
    if inCache('getGeometry@'+sparams):
        return fromCache('getGeometry@'+sparams)
    errorMsg, params = validateClippingResolution(sparams)
    if (errorMsg != ''):
        raise Exception('SISMOI Err: getGeometry: ' + errorMsg + '\nParams: '+sparams)
    if not params['resolution'] in mapTemplates:
        cur.execute("SELECT geojson FROM geojson WHERE name = '{0}'".format(params['resolution']))
        mapTemplates[params['resolution']]=cur.fetchone()[0]
        map=json.loads(mapTemplates[params['resolution']])
        try:
            i=0
            while i < len(map['features']):
                if params['clipping'] != map['features'][i]['properties']['state']:
                    del map['features'][i]
                else:
                    i+=1
        except Exception as e:
            print(e)
        ret=json.dumps(map)
    else:
        ret=mapTemplates[params['resolution']]
    toCache('getGeometry@'+sparams,ret)
    return ret

@app.route('/sismoi/getMapData/<sparams>', methods=['GET'])
def getMapData(sparams):
    if inCache('getMapData@'+sparams):
        return fromCache('getMapData@'+sparams)
    errorMsg, params = validateParams(sparams)
    if (errorMsg != ''):
        raise Exception('getMapData: ' + errorMsg)
    data = None
    if params['resolution'] == 'municipio':
        data = getIndicatorByCounty(params)
    elif params['resolution'] == 'microrregiao':
        data = getIndicatorByMicroregion(params)
    elif params['resolution'] == 'mesorregiao':
        data = getIndicatorByMesoregion(params)
    elif params['resolution'] == 'estado':
        data = getIndicatorByState(params)
    if len(data) == 0:
        raise Exception('SISMOI Err: getMapData: Não existem dados para esses parâmetros: {0}'.format(sparams))
    ret=json.dumps(data)
    toCache('getMapData@'+sparams,ret)
    return ret

@app.route('/sismoi/getGeometryAndData/<sparams>', methods=['GET'])
def getGeometryAndData(sparams):
    if inCache('getGeometryAndData@'+sparams):
        return fromCache('getGeometryAndData@'+sparams)
    errorMsg, params = validateParams(sparams)
    if (errorMsg != ''):
        raise Exception('getGeometryAndData: ' + errorMsg)
    geometry = json.loads(getGeometry(sparams))
    mapdata = json.loads(getMapData(sparams))
    keyfield = 'state' if params['resolution'] == 'estado' else 'id'
    for feature in geometry['features']:
        row = findElement(mapdata, keyfield, feature)
        feature['properties']['value'] = float(row['value'])
        feature['properties']['style'] = {'color': '#d7d7d7',
                                          'fillcolor': featureColor(float(row['value']), row['pessimist']), 'weight': 1}
    ret=json.dumps(geometry)
    toCache('getGeometryAndData@'+sparams,ret)
    return ret

@app.route('/sismoi/getTotal/<sparams>', methods=['GET'])
def getTotal(sparams):
    if inCache('getTotal@'+sparams):
        return fromCache('getTotal@'+sparams)
    errorMsg, params = validateParams(sparams)
    if (errorMsg != ''):
        raise Exception('getTotal: ' + errorMsg)
    if params['resolution'] == 'municipio':
        data = getTotalByCounty(params)
    elif params['resolution'] == 'microrregiao':
        data = getTotalByMicroregion(params)
    elif params['resolution'] == 'mesorregiao':
        data = getTotalByMesoregion(params)
    elif params['resolution'] == 'estado':
        data = getTotalByState(params)
    else:
        raise Exception('SISMOI Err: Parâmetro resolution inválido: {0}',params['resolution'])
    ret=json.dumps(data)
    toCache('getTotal@'+sparams,ret)
    return ret

if __name__ == "__main__":
    try:
        conn, cur = connect()
        clippings = clippings + getStates()
        years = getYears()
        args = ProcessCmdLine()
        cacheType=args.cachetype
        start=datetime.datetime.now()

        app.run(host=args.host, port=args.port, debug=args.debug)
    except Exception as e:
        print(e)
        exit(-1)
