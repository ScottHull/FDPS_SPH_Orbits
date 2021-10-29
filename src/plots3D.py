cube_verts = [
    [
        [0, 0, 0], [0, 0, 1]
    ],
[
        [0, 0, 0], [1, 0, 0]
    ],
[
        [1, 0, 0], [1, 0, 1]
    ],
[
        [0, 0, 0], [0, 1, 0]
    ],
[
        [0, 0, 1], [0, 1, 1]
    ],
[
        [0, 0, 1], [1, 0, 1]
    ],
[
        [1, 0, 1], [1, 1, 1]
    ],
[
        [1, 1, 0], [1, 1, 1]
    ],
[
        [1, 0, 0], [1, 1, 0]
    ],
[
        [0, 1, 0], [0, 1, 1]
    ],
[
        [0, 1, 1], [1, 1, 1]
    ],
[
        [0, 1, 0], [1, 1, 0]
    ],
]

def get_cube_verts(square_scale):
    verts = []
    for i in cube_verts:
        tmp = []
        for j in i:
            tmp.append([(-square_scale if k == 0 else square_scale) for k in j])
        verts.append(list(zip(tmp[0], tmp[1])))
    return verts