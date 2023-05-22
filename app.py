import json
from flask import Flask, jsonify, request
from flask_restful import Api, Resource
from flask_cors import CORS
#import scipy as optimize
#import scipy.integrate as quad
# import pyrebase
# import urllib3
# from firebase_admin import credentials, db
# import firebase_admin
import os
import numpy as np
import matplotlib.pyplot as plt

app = Flask(__name__)
cors = CORS(app, resources={r"/v1/api/*": {"origins": "*"}})
api = Api(app)
# cred = credentials.Certificate('key.json')
# firebase_admin.initialize_app(cred, {
#     'databaseURL': 'https://watermeasurement-8bdae-default-rtdb.firebaseio.com/'
# })
# firebase = pyrebase.initialize_app(confirg)
# db = firebase.database()

# class Database(Resource):
#     def get(self):
#         Users = db.reference('users')
#         check = Users.child('user3')
#         print(check.get())
#         dict = {
#             "email": 'eashamondal01@gmail.com',
#             'password': 'Easha@2001'
#         }
#         email = 'eashamondal01@gmail.com'
#         dict2 = {}
#         dict2[email] = dict
#         print(dict2)
#         Users = Users.child()
#         Users.set(dict)
#         # Users.push(dict2)

# api.add_resource(Database, '/v1/api/datas')

def error(message):
    result = {
        "error": message,
        "status": False,
    }
    return jsonify(result)

# class Login(Resource):
#     def post(self):
#         data = request.get_json()
#         name = data["name"]
#         password = data["password"]
#         dict = Users.get().val()
#         json_string = json.dumps(dict)
#         obj = json.loads(json_string)
#         user_list = obj.values()
#         for user in user_list:
#             if user['name'] == name and user['password'] == password:
#                 result = {
#                     "error": "",
#                     "status": True,
#                     "name": name,
#                 }
#                 return jsonify(result)
#         return error("Wrong username or password")


def Activated_sludge(time, volume_of_reactor, recycle_ratio, area_of_clarifier, height_of_clarifier, bod_cod_ratio, BOD, iCOD, mlvss, flowrate, index):
    def RK4(func,yi,ti,dt):
        k1 = dt*func(ti,yi)
        k2 = dt*func(ti+0.5*dt,yi+0.5*k1)
        k3 = dt*func(ti+0.5*dt,yi+0.5*k2)
        k4 = dt*func(ti+dt,yi+k3)
        return yi+(1./6*(k1+2*k2+2*k3+k4))
    #defining time dependent differential   
    def dConc(t,C):
        dC = np.zeros(6) #enter no.of differential equations
        #write all differential equation using parameters as variables
        #mass balance on soluble inerts
        dC[0] = Q1*Cin[0] - Q0*C[0]
        #mass balance of readily biodegradable
        dC[1] = Q1*Cin[1] - Q0*C[1] - ((1/YH)*muh*(C[1]/(Ks + C[1]))*C[4]) + (kh*((C[3]/C[4])/(Kx+(C[3]/C[4])))*C[4])
        #mass balance on insoluble inerts
        dC[2] = Q1*Cin[2] - Q0*C[2]
        # mass balance on slowly biodegradable
        dC[3] = Q1*Cin[3] - Q0*C[3] + ((1-fp)*(bh*C[4])) - (kh*((C[3]/C[4])/(Kx+(C[3]/C[4])))*C[4])
        #mass balance on active sludge
        dC[4] = Q1*Cin[4] - Q0*C[4] + ((1/YH)*muh*(C[1]/(Ks+C[1]))*C[4]) - (bh*C[4])
        #mass balance on inactive sludge
        dC[5] = Q1*Cin[5] - Q0*C[5] + (fp*bh*C[4])
        return dC
    #creating time steps
    t=144 
    tmax=time/24
    dT=0.001
    nT = (tmax/dT)+1
    nT = int(nT)
    T = np.linspace(0,tmax,nT)
    #entering the value of parameters
    rr=recycle_ratio
    BOD=BOD 
    COD=iCOD 
    YH, muh, kh,Ks, fp, Kx, bh = 0.6, 0.9, 2.5, 8, 0.4, 0.7, 0.82
    Qf=flowrate 
    Qr=rr*Qf
    Volume = volume_of_reactor 
    MLVSS = mlvss 
    Q0=(Qr+Qf)/Volume
    Q1=Q0
    #Q0,Q1=0.3,0.3
    #entering inlet concn
    Cin = [(0.5*(COD-BOD)),(0.8*BOD),(0.2*BOD),(0.5*(COD-BOD)),0,0]
    #Cin=np.array(Cin)
    #Cin = [170, 420, 35, 1675, 0, 0]
    T0 = 0.0
    Concn = np.zeros((nT,6))
    #entering initial concn
    #Cinit = [250, 650, 60, 2500, 100, 0]
    Cinit = [(0.5*(COD-BOD)),(0.8*BOD),(0.2*BOD),(0.5*(BOD-COD)),(0.95*MLVSS),(0.05*MLVSS)]
    print(Cinit)
    #Cinit=np.array(Cinit)
    Concn[0,:] = Cinit
    BOD=np.zeros((nT))
    #computing values of components at different time steps
    #calculating oxygen required
    OR = np.zeros((nT,1))
    for iT in range(1,nT):
        Concn[iT,:]= RK4(dConc,Concn[iT-1,:],T[iT-1],dT)
        OR[iT] = ((1-YH)/YH)*muh*(Concn[iT,1]/(Ks+Concn[iT,1])*Concn[iT,4])
    BOD=Concn[:,1]+Concn[:,3]
    COD=Concn[:,0]+Concn[:,2]+Concn[:,1]+Concn[:,3]
    MLVSS=Concn[:,5]+Concn[:,4]
    TC=8*10**4*MLVSS
    TSS=Concn[:,2]+Concn[:,3]+Concn[:,4]+Concn[:,5]
    BOD = BOD[::600].tolist()
    COD = COD[::600].tolist()
    TSS = TSS[::600].tolist()
    TC = TC[::600].tolist()
    result = {
        'bod': int(BOD[len(BOD)-1]),
        'cod': int(COD[len(COD)-1]),
        'tss': int(TSS[len(TSS)-1]),
        'total_caliform': int(TC[len(TC)-1]),
        '_BOD': BOD,
        '_label_BOD': list(range(1, len(BOD)+1)),
        '_COD': COD,
        '_label_COD': list(range(1, len(COD)+1)),
        '_TSS': TSS,
        '_label_TSS': list(range(1, len(TSS)+1)),
        '_TC': TC,
        '_label_TC': list(range(1, len(TC)+1)),
        'time': time
    }
    return result
# def Activated_sludge(time, volume_of_reactor, recycle_ratio, area_of_clarifier, height_of_clarifier, bod_cod_ratio, BOD, iCOD, mlvss, flowrate, index):
#     YH, muh, Ks, fp, Kx, bh, kh, Q1, Q0 = 0.6, 0.8, 8, 0.4, 0.7, 0.82, 4.5, 0.312, 0.312
#     tmax = time/24
#     dT = 0.01
#     nT = (tmax/dT) + 1
#     nT = int(nT)
#     T = np.linspace(0, tmax, nT)
#     Qf = flowrate
#     Reratio = recycle_ratio
#     Volume = volume_of_reactor
#     Qr = Reratio*Qf
#     Q1 = Qr + Qf
#     Q1 = (Volume/Q1)
#     Q0 = Q1
#     T0 = 0
#     def RK4(func, yi, ti, dt):
#         k1 = dt*func(ti, yi)
#         k2 = dt*func(ti+0.5*dt, yi+0.5*k1)
#         k3 = dt*func(ti+0.5*dt, yi+0.5*k2)
#         k4 = dt*func(ti+dt, yi+k3)
#         kil = (1./6*(k1+2*k2+2*k3+k4))
#         return yi+(1./6*(k1+2*k2+2*k3+k4))
#     def dConc(t, C):
#         dC = np.zeros(6)
#         Q0 = Q1
#         Cin = [
#             0.5*(BOD-iCOD),
#             0.8*BOD,
#             0.2*BOD,
#             0.5*(iCOD-BOD),
#             0, 0
#         ]
#         # Cin[0] = 0.5*(COD - BOD)
#         # Cin[1] = 0.8*BOD
#         # Cin[2] = 0.2*BOD
#         # Cin[3] = 0.5*(COD - BOD)
#         # Cin[4] = 0.95*MLVSS Cinit, Cin=0
#         # Cin[5] = 0.05*MLVSS Cinit Cin = 0
#         dC[0] = Q1*Cin[0] - Q0*C[0]
#         dC[1] = Q1*Cin[1] - Q0*C[1] - ((1/YH)*muh*(C[1]/(Ks + C[1]))*C[4]) + (kh*((C[3]/C[4])/(Kx+(C[3]/C[4])))*C[4])
#         dC[2] = Q1*Cin[2] - Q0*C[2]
#         dC[3] = Q1*Cin[3] - Q0*C[3] + ((1 - fp)*(bh*C[4])) - (kh*((C[3]/C[4])/(Kx+(C[3]/C[4])))*C[4])
#         dC[4] = Q1*Cin[4] - Q0*C[4] + ((1/YH)*muh*(C[1]/(Ks + C[1]))*C[4]) - (bh*C[4])
#         dC[5] = Q1*Cin[5] - Q0*C[5] + (fp*bh*C[4])
#         return dC
#     Concn = np.zeros((nT, 6))
#     print(type(Concn))
#     Cinit = [
#         0.5*(iCOD-BOD),
#         0.8*BOD,
#         0.2*BOD,
#         0.5*(iCOD-BOD),
#         0.95*mlvss,
#         0.05*mlvss
#     ]
#     Concn[0, :] = Cinit
#     UltimateBOD = np.zeros((nT, 1))
#     COD = np.zeros((nT, 1))
#     ReadilyB = np.zeros((nT, 1))
#     slowlyB = np.zeros((nT, 1))
#     SolInerts = np.zeros((nT, 1))
#     Insolinerts = np.zeros((nT, 1))
#     Activatesludges = np.zeros((nT, 1))
#     Inactivesludge = np.zeros((nT, 1))
#     TSS = np.zeros((nT, 1))
#     TC = np.zeros((nT, 1))
#     OR = np.zeros((nT, 1))
#     for iT in range(1, nT):
#         Concn[iT, :] = RK4(dConc, Concn[iT-1, :], T[iT-1], dT)
#         # OR[iT] = ((1-YH)/YH)*muh*(Concn[iT, 1]/(Ks+Concn[iT, 1])*Concn[iT, 4])
#         ReadilyB[iT] = Concn[iT, 1]
#         slowlyB[iT] = Concn[iT,3]
#         SolInerts[iT] = Concn[iT,0]
#         Insolinerts[iT] = Concn[iT,2]
#         Activatesludges[iT] = Concn[iT,4]
#         Inactivesludge[iT] = Concn[iT,5]
#     UltimateBOD[iT] = Concn[iT-1,1]+Concn[iT-1, 3]
#     COD[iT] = Concn[iT,0]+Concn[iT,1]+Concn[iT-1,2]+Concn[iT-1,3]
#     TSS[iT] = Concn[iT,3]+Concn[iT,2]+Concn[iT,4]+Concn[iT,5]
#     TC[iT] = Concn[iT,4]
#     index = int(index)
#     print(len(UltimateBOD))
#     # bod_label = graph_for_bod
#     graph_for_bod = UltimateBOD[::50].tolist()
#     print(len(graph_for_bod))
#     graph_for_bod = [int(i[0]) for i in graph_for_bod]
#     graph_for_bod[0] = BOD
#     graph_for_cod = COD[::50].tolist()
#     print(len(graph_for_cod))
#     graph_for_cod = [int(i[0]) for i in graph_for_cod]
#     graph_for_cod[0] = iCOD
#     graph_for_TSS = TSS[::50].tolist()
#     graph_for_TSS = [int(i[0]) for i in graph_for_TSS]
#     graph_for_TSS[0] = 100000
#     graph_for_TC = TC[::50].tolist()
#     graph_for_TC = [int(i[0]) for i in graph_for_TC]
#     graph_for_TC[0] = 4000
#     result = {
#         'bod': int(UltimateBOD[len(UltimateBOD)-1][0]),
#         'cod': int(COD[len(COD)-1][0]),
#         'tss': int(TSS[len(TSS)-1][0]),
#         'total_caliform': int(TC[len(TC)-1][0]),
#         '_BOD': graph_for_bod,
#         '_label_BOD': list(range(1, len(UltimateBOD)+1)[::50]),
#         '_COD': graph_for_cod,
#         '_label_COD': list(range(1, len(graph_for_cod)+1)),
#         '_TSS': graph_for_TSS,
#         '_label_TSS': list(range(1, len(graph_for_TSS)+1)),
#         '_TC': graph_for_TC,
#         '_label_TC': list(range(1, len(graph_for_TC)+1)),
#         'time': time
#     }
#     print(type(UltimateBOD[len(UltimateBOD)-1][0]))
#     return result

def Hydrodynamic_cavitation(bod, cod, batch_volume, number_of_passes, time_cavitation, pH, pressure_cavitation):
    total = batch_volume + number_of_passes + pH + pressure_cavitation
    tss = total*time_cavitation
    total_caliform = total + pH*pressure_cavitation
    result = {
        'bod': total*time_cavitation + bod,
        'cod': total*time_cavitation + cod,
        'tss': tss,
        'total_caliform': total_caliform,
        'time': time_cavitation
    }
    return result

def Disinfection(tc, time_disinfection, concentration_disinfection, bod, cod, tss):
    n, m, k1, k = 0.22, 0.2, 0.08, 1.2
    def func(x):
        return ((np.exp(-k1*n*x))*(x**(m-1)))
    def func1(x):
        return ((np.exp(-k1*n*x)))
    C0 = concentration_disinfection
    # N0 = 3*1000000
    N0 = tc
    ti = 0
    tmax = time_disinfection
    nT = 1000
    tf = np.linspace(0.001, tmax, nT)
    val = np.linspace(0.001, tmax, nT)
    CCl = np.linspace(0.001, tmax, nT)
    CoutCl = np.linspace(0.001, tmax, nT)
    for iT in range(0, nT):
        estimate = quad.quadrature(func, ti, tf[iT], maxiter=100, tol=0.001)
        val[iT] = estimate[0]
        estimate1 = quad.quadrature(func1, ti, tf[iT], maxiter=100, tol=0.001)
        CCl[iT] = C0*estimate1[0]
        CoutCl[iT] = C0*np.exp(-k1*tf[iT])
    NCl = N0*(10**(-k*m*(C0**n)*val))
    NCl = NCl[::78].tolist()
    NCl = [int(x) for x in NCl]
    print(len(NCl), 'hi')
    # logCl = np.log(N0/NCl)
    result = {
        '_TC': NCl,
        '_label_TC': list(range(1, len(NCl)+1)),
        'time': time_disinfection
    }
    return result




# class SignUp(Resource):
#     def post(self):
#         data = request.get_json()
#         name = data["name"]
#         password = data["password"]
#         dict = db.child('users').get().val()
#         json_string = json.dumps(dict)
#         obj = json.loads(json_string)
#         user_list = obj.values()
#         for user in user_list:
#             if user['name'] == name:
#                 return error('User already exists with this username')
#         new_user = {
#             "name": name,
#             "password": password
#         }
#         db.child('users').push(new_user)
#         result = {
#             "error": "",
#             "status": True,
#             "name": name
#         }
#         return jsonify(result)



class Test(Resource):
    def get(self):
        result = {
            "message": "This is a test api",
        }
        return jsonify(result)

class OutputEstimations(Resource):
    def post(self, index):
        data = request.get_json()
        if "flowrate" in data.keys():
            flowrate = int(data["flowrate"])
        else:
            return error("Flowrate is required")
        if "res_time_sludge" in data.keys():
            res_time_sludge = int(data["res_time_sludge"])
        else:
            return error("Time for Activated_sludge is required")
        if "volume_of_reactor" in data.keys():
            volume_of_reactor = int(data["volume_of_reactor"])
        else:
            return error("Volume of Reactor is required")
        if "recycle_ratio" in data.keys():
            recycle_ratio = float(data["recycle_ratio"])
        else:
            return error("Recycle ratio is required")
        area_of_clarifier = int(data["area_of_clarifier"])
        height_of_clarifier = int(data["height_of_clarifier"])
        bod_cod_ratio = int(data["bod_cod_ratio"])
        if "bod" in data.keys():
            bod = int(data["bod"])
        else:
            return error("Bod is required")
        if "cod" in data.keys():
            cod = int(data["cod"])
        else:
            return error("Cod is required")
        mlvss = int(data["mlvss"])
        batch_volume = 0
        number_of_passes = 0
        time_cavitation = 0
        pH = 0
        pressure_cavitation = 0
        if "time_disinfection" in data.keys():
            time_disinfection = int(data["time_disinfection"])
        else:
            return error("TimeDisinfection is required")
        if "concentration_disinfection" in data.keys():
            concentration_disinfection = int(data["concentration_disinfection"])
        else:
            return error("ConcentrationDisinfection is required")
        result_activated_sludge = Activated_sludge(res_time_sludge, volume_of_reactor, recycle_ratio, area_of_clarifier, height_of_clarifier, bod_cod_ratio, bod, cod, mlvss, flowrate, index)
        result_hydrodynamic_cavitation = Hydrodynamic_cavitation(result_activated_sludge['bod'], result_activated_sludge['cod'], batch_volume, number_of_passes, time_cavitation, pH, pressure_cavitation)
        result_disinfection = Disinfection(result_activated_sludge['total_caliform'], time_disinfection, concentration_disinfection, bod, cod, result_activated_sludge['tss'])
        list = {
            'result_activated_sludge': result_activated_sludge,
            'result_hydrodynamic_cavitation': result_hydrodynamic_cavitation,
            'result_disinfection': result_disinfection,
        }
        result = {
            "error": "",
            "status": True,
            "result": list,
            "total_time": res_time_sludge+time_cavitation+time_disinfection
        }
        return jsonify(result)

class Pre_Calculation(Resource):
    def post(self):
        data = request.get_json()
        bod = int(data['bod'])
        cod = int(data['cod'])
        time = int(data['time'])
        Tss = int(data['tss'])
        tc = int(data['tc'])
        t = range(1, time)[::time//10]
        TSS = []
        BOD = [0]
        COD = [0]
        TC = [0]
        for i in t:
            TSS.append(int((1-(0.4*i/time))*Tss))
            BOD.append(bod)
            COD.append(cod)
            TC.append(tc)
        result = {
            'BOD': BOD,
            'COD': COD,
            'TC': TC,
            'TSS': TSS,
            'label': list(t)
        }
        return jsonify(result)

class Terticary_calculation(Resource):
    def post(self):
        data = request.get_json()
        bod = int(data['bod'])
        cod = int(data['cod'])
        time = int(data['time'])
        tss = int(data['tss'])
        tc = int(data['tc'])
        t = range(1, time)[::time//10]
        TSS = [0]
        BOD = []
        COD = []
        TC = [0]
        for i in t:
            BOD.append(int((1-(0.99*i/time))*bod))
            COD.append(int((1-(0.85*i/time))*cod))
            TSS.append(tss)
            TC.append(tc)
        t = list(t)
        result = {
            'BOD': BOD,
            'COD': COD,
            'TC': TC,
            'TSS': TSS,
            'label': t
        }
        return jsonify(result)

class Disinfection_calculation(Resource):
    def post(self):
        data = request.get_json()
        bod = int(data['bod'])
        cod = int(data['cod'])
        tss = int(data['tss'])
        tc = int(data['tc'])
        time = int(data['time_disinfection'])
        conc = int(data['concentration_disinfection'])
        dict = Disinfection(tc, time, conc, bod, cod, tss)
        BOD = [0]
        COD = [0]
        TSS = [0]
        TC = []
        label = dict['_label_TC']
        for i in label:
            BOD.append(bod)
            COD.append(cod)
            TSS.append(tss)
        result = {
            'BOD': BOD,
            'COD': COD,
            'TC': dict['_TC'],
            'TSS': TSS,
            'label': dict['_label_TC']
        }
        return jsonify(result)

api.add_resource(Disinfection_calculation, '/v1/api/Disinfection_calculation')
api.add_resource(Terticary_calculation, '/v1/api/Terticary_calculation')
api.add_resource(Pre_Calculation, '/v1/api/Pre_Calculation')
api.add_resource(Test, '/v1/api/test')
api.add_resource(OutputEstimations, '/v1/api/outputEstimations/<index>')

if __name__ == '__main__':
    port = int(os.environ.get("PORT", 5000))
    app.run(debug=True, port= port, host='0.0.0.0')
