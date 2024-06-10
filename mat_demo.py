#!/usr/bin/env python

import Geo
# Import the matlab module only after you have imported
# MATLAB Compiler SDK generated Python modules.
import matlab
import redis
import json

# Connect to Redis server
r = redis.Redis(host='localhost', port=6379, db=0)

# Add items to the queue
queue_name = 'matlab_data'

# filename = 'data/bent_pipe_closed_lr.off'
filename = 'data/spot_rr.off'

def load_off_file(filename):
    with open(filename, 'r') as file:
        if 'OFF' != file.readline().strip():
            raise ValueError('Not a valid OFF header')

        n_verts, n_faces, _ = map(int, file.readline().strip().split(' '))

        vertices = []
        for _ in range(n_verts):
            vertices.extend(map(float, file.readline().strip().split(' ')))

        vertices = [vertices[i:i + 3] for i in range(0, len(vertices), 3)]

        faces = []
        for _ in range(n_faces):
            face = list(map(int, file.readline().strip().split(' ')))
            faces.append(face[1:])

    return vertices, faces

try:
    my_Geo = Geo.initialize()
except Exception as e:
    print('Error initializing Geo package\\n:{}'.format(e))
    exit(1)

try:
    vertices, faces = load_off_file(filename)

    vIn = matlab.double(vertices)
    # Adjust indices from 0-based to 1-based
    faces_matlab = [[index + 1 for index in face] for face in faces]
    fIn = matlab.double(faces_matlab)

    x0In = matlab.double([2273], size=(1, 1))
    given_vf_facesIn = matlab.double([4736, 2703], size=(1, 2))
    given_vf_valsIn = matlab.double([1.6256, 1.6952, -0.3518, 0.3193, -0.6234, 0.0335], size=(2, 3))
    u_vfaOut = my_Geo.geo(vIn, fIn, x0In, given_vf_facesIn, given_vf_valsIn)
    print(u_vfaOut, sep='\\n')

    print("Dimensions of u_vfaOut:", u_vfaOut.size)
    print("val: ", u_vfaOut[0])
    gd = []
    for elem in u_vfaOut:
        gd.append(elem[0])
    print("Dimensions of gd:", len(gd))

    data = {
        'vertices': vertices,
        'faces': faces,
        'u_vfaOut': gd
    }
    data_json = json.dumps(data)
    r.lpush(queue_name, data_json)

except Exception as e:
    print('Error occurred during program execution\\n:{}'.format(e))

my_Geo.terminate()
