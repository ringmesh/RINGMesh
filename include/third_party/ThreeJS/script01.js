function init() {
    scene = new THREE.Scene();
    renderer = new THREE.WebGLRenderer({
        alpha: false
    });
    renderer.setClearColor(0xffffff, 1.0);
    renderer.setSize(window.innerWidth, window.innerHeight);

    document.body.appendChild(renderer.domElement);
    camera = new THREE.PerspectiveCamera(40, window.innerWidth / window.innerHeight, 0.1, 100000);
    initControl();
    initMaterials();
    initLights();
    loadObjects();
    setCameraPlace();
    scene.add(camera)
};

function animate() {
    requestAnimationFrame(animate);
    render();
    update()
};

function update() {
    controls.update();
    stats.update()
};

function render() {
    light.position.set(camera.position.x, camera.position.y, camera.position.z);
    renderer.render(scene, camera)
};

function setCameraPlace() {
    if (meshes.length == 0) {
        return
    }
    meshes[0].geometry.computeBoundingBox();
    minX = meshes[0].geometry.boundingBox.min.x;
    maxX = meshes[0].geometry.boundingBox.max.x;
    minY = meshes[0].geometry.boundingBox.min.y;
    maxY = meshes[0].geometry.boundingBox.max.y;
    minZ = meshes[0].geometry.boundingBox.min.z;
    maxZ = meshes[0].geometry.boundingBox.max.z;
    for (obj = 0; obj < meshes.length; ++obj) {
        if (meshes[obj].visible) {
            mesh = meshes[obj];
            mesh.geometry.computeBoundingBox();
            box3 = mesh.geometry.boundingBox;
            minX = Math.min(minX, box3.min.x);
            minY = Math.min(minY, box3.min.y);
            minZ = Math.min(minZ, box3.min.z);
            maxX = Math.max(maxX, box3.max.x);
            maxY = Math.max(maxY, box3.max.y);
            maxZ = Math.max(maxZ, box3.max.z)
        }
    }



    center = new THREE.Vector3((minX + maxX) * 0.5, (minY + maxY) * 0.5, (minZ + maxZ) * 0.5);
    distance = Math.sqrt((maxX - minX) * (maxX - minX) + (maxY - minY) * (maxY - minY) + (maxZ - minZ) * (maxZ - minZ));
    for (obj = 0; obj < meshes.length; ++obj) {
        if (meshes[obj].visible) {
            meshes[obj].translateX(-center.x);
            meshes[obj].translateY(-center.y);
            meshes[obj].translateZ(-center.z);
        }
    }
    camera.position.set(distance * 0.75, distance * 0.75, -distance * 0.75);
    camera.up = new THREE.Vector3(0, 0.0, -1.0)
};

function initControl() {
    controls = new THREE.TrackballControls(camera, renderer.domElement);
    controls.addEventListener('change');
    controls.rotateSpeed *= 2.0;
    controls.zoomSpeed *= 2.0;
    controls.panSpeed *= 2.0
};

function initMaterials() {
    material = new THREE.MeshLambertMaterial({
        color: 0xfffff,
        shading: THREE.SmoothShading,
        side: THREE.DoubleSide
    })
};

function initLights() {
    ambientlight = new THREE.AmbientLight(0x404040);
    scene.add(ambientlight);
    light = new THREE.PointLight(0xfffaff, 0.8, 0);
    light.position.set(camera.position.x, camera.position.y, camera.position.z);
    scene.add(light)
};



function loadTSurf(numbers) {
    number_of_points = numbers[0];
    number_of_triangles = numbers[1];
    var geometry = new THREE.Geometry();
    for (var p = 0; p < number_of_points; p++) {
        var v = new THREE.Vector3(numbers[2 + p * 3 + 0], numbers[2 + p * 3 + 1], numbers[2 + p * 3 + 2]);
        geometry.vertices.push(v)
    }
    var offset = 3 * number_of_points + 2;
    for (var t = 0; t < number_of_triangles; t++) {
        var f = new THREE.Face3(numbers[offset + t * 3 + 0] , numbers[offset + t * 3 + 1] , numbers[offset + t * 3 + 2] );
        geometry.faces.push(f);
    }

    geometry.computeFaceNormals();
    mat = new THREE.MeshLambertMaterial(material);
    var surface = new THREE.Mesh(geometry, mat);
    meshes.push(surface);

};

function loadPLine(numbers){

	var type = THREE.LineStrip ; //THREE.LinePieces
	var geometry = new THREE.Geometry();
	for (var p = 0; p < numbers.length / 3 ; p++) {
		var point = new THREE.Vector3(numbers[3*p+0],numbers[3*p+1],numbers[3*p+2]);
		geometry.vertices.push(point);
	}
	var line = new THREE.Line( geometry, new THREE.LineBasicMaterial( { color: Math.random() * 0xffffff , linewidth:5 }  ), type );
	 meshes.push(line);
}

function loadTSolid(numbers){
	loadTSurf(numbers);
}

function loadVSet(numbers){

	var geometry = new THREE.Geometry();
	for (var p = 0; p < numbers.length / 3 ; p++) {
		var point = new THREE.Vector3(numbers[3*p+0],numbers[3*p+1],numbers[3*p+2]);
		geometry.vertices.push(point);
	}
	var pointCloud = new THREE.PointCloud( geometry, new THREE.PointCloudMaterial( { color: Math.random() * 0xffffff  } ) );
	 meshes.push(pointCloud);

}


function addMeshes(){
            for (m = 0; m < meshes.length; ++m) {
                scene.add(meshes[m]);
            }
}
