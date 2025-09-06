//#include "C:/Wichtig/System/Static/Library/WindowEngine1.0.h"
#include "/home/codeleaded/System/Static/Library/WindowEngine1.0.h"
#include "./Math3D.h"

mesh meshCube;
mat4x4 matProj;
vec3d vCamera;
vec3d vLookDir;
float fYaw;
float fTheta;

Sprite sprTex1;
float *pDepthBuffer;

void MakeCube(vec3d p,vec3d d,Pixel c){
	triangle tris[12] = {
		// SOUTH
		{ 0.0f, 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 0.0f, 1.0f,    1.0f, 1.0f, 0.0f, 1.0f,		0.0f, 1.0f, 1.0f,		0.0f, 0.0f, 1.0f,		1.0f, 0.0f, 1.0f, c}, 
		{ 0.0f, 0.0f, 0.0f, 1.0f,    1.0f, 1.0f, 0.0f, 1.0f,    1.0f, 0.0f, 0.0f, 1.0f,		0.0f, 1.0f, 1.0f,		1.0f, 0.0f, 1.0f,		1.0f, 1.0f, 1.0f, c},
						  																			   
		// EAST           																			   
		{ 1.0f, 0.0f, 0.0f, 1.0f,    1.0f, 1.0f, 0.0f, 1.0f,    1.0f, 1.0f, 1.0f, 1.0f,		0.0f, 1.0f, 1.0f,		0.0f, 0.0f, 1.0f,		1.0f, 0.0f, 1.0f, c},
		{ 1.0f, 0.0f, 0.0f, 1.0f,    1.0f, 1.0f, 1.0f, 1.0f,    1.0f, 0.0f, 1.0f, 1.0f,		0.0f, 1.0f, 1.0f,		1.0f, 0.0f, 1.0f,		1.0f, 1.0f, 1.0f, c},
						   																			   
		// NORTH           																			   
		{ 1.0f, 0.0f, 1.0f, 1.0f,    1.0f, 1.0f, 1.0f, 1.0f,    0.0f, 1.0f, 1.0f, 1.0f,		0.0f, 1.0f, 1.0f,		0.0f, 0.0f, 1.0f,		1.0f, 0.0f, 1.0f, c},
		{ 1.0f, 0.0f, 1.0f, 1.0f,    0.0f, 1.0f, 1.0f, 1.0f,    0.0f, 0.0f, 1.0f, 1.0f,		0.0f, 1.0f, 1.0f,		1.0f, 0.0f, 1.0f,		1.0f, 1.0f, 1.0f, c},
						   																			   
		// WEST            																			   
		{ 0.0f, 0.0f, 1.0f, 1.0f,    0.0f, 1.0f, 1.0f, 1.0f,    0.0f, 1.0f, 0.0f, 1.0f,		0.0f, 1.0f, 1.0f,		0.0f, 0.0f, 1.0f,		1.0f, 0.0f, 1.0f, c},
		{ 0.0f, 0.0f, 1.0f, 1.0f,    0.0f, 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 0.0f, 1.0f,		0.0f, 1.0f, 1.0f,		1.0f, 0.0f, 1.0f,		1.0f, 1.0f, 1.0f, c},
						   																			   
		// TOP             																			   
		{ 0.0f, 1.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f, 1.0f,    1.0f, 1.0f, 1.0f, 1.0f,		0.0f, 1.0f, 1.0f,		0.0f, 0.0f, 1.0f,		1.0f, 0.0f, 1.0f, c},
		{ 0.0f, 1.0f, 0.0f, 1.0f,    1.0f, 1.0f, 1.0f, 1.0f,    1.0f, 1.0f, 0.0f, 1.0f,		0.0f, 1.0f, 1.0f,		1.0f, 0.0f, 1.0f,		1.0f, 1.0f, 1.0f, c},
						   																			  
		// BOTTOM          																			  
		{ 1.0f, 0.0f, 1.0f, 1.0f,    0.0f, 0.0f, 1.0f, 1.0f,    0.0f, 0.0f, 0.0f, 1.0f,		0.0f, 1.0f, 1.0f,		0.0f, 0.0f, 1.0f,		1.0f, 0.0f, 1.0f, c},
		{ 1.0f, 0.0f, 1.0f, 1.0f,    0.0f, 0.0f, 0.0f, 1.0f,    1.0f, 0.0f, 0.0f, 1.0f,		0.0f, 1.0f, 1.0f,		1.0f, 0.0f, 1.0f,		1.0f, 1.0f, 1.0f, c},
	};

	for(int i = 0;i<12;i++){
		for(int j = 0;j<3;j++){
			tris[i].p[j] = vec3d_Add(p,vec3d_Make(tris[i].p[j].x * d.x,tris[i].p[j].y * d.y,tris[i].p[j].z * d.z));
		}
		Vector_Push(&meshCube.tris,&tris[i]);
	}
}

void MakePlane(vec3d p,vec3d d,int Plane,Pixel c){
	triangle tris[12] = {
		// SOUTH
		{ 0.0f, 0.0f, 0.0f, 1.0f,    0.0f, 1.0f, 0.0f, 1.0f,    1.0f, 1.0f, 0.0f, 1.0f,		0.0f, 1.0f, 1.0f,		0.0f, 0.0f, 1.0f,		1.0f, 0.0f, 1.0f, c}, 
		{ 0.0f, 0.0f, 0.0f, 1.0f,    1.0f, 1.0f, 0.0f, 1.0f,    1.0f, 0.0f, 0.0f, 1.0f,		0.0f, 1.0f, 1.0f,		1.0f, 0.0f, 1.0f,		1.0f, 1.0f, 1.0f, c},
						  																			   
		// EAST           																			   
		{ 1.0f, 0.0f, 0.0f, 1.0f,    1.0f, 1.0f, 0.0f, 1.0f,    1.0f, 1.0f, 1.0f, 1.0f,		0.0f, 1.0f, 1.0f,		0.0f, 0.0f, 1.0f,		1.0f, 0.0f, 1.0f, c},
		{ 1.0f, 0.0f, 0.0f, 1.0f,    1.0f, 1.0f, 1.0f, 1.0f,    1.0f, 0.0f, 1.0f, 1.0f,		0.0f, 1.0f, 1.0f,		1.0f, 0.0f, 1.0f,		1.0f, 1.0f, 1.0f, c},
						   																			   
		// NORTH           																			   
		{ 1.0f, 0.0f, 1.0f, 1.0f,    1.0f, 1.0f, 1.0f, 1.0f,    0.0f, 1.0f, 1.0f, 1.0f,		0.0f, 1.0f, 1.0f,		0.0f, 0.0f, 1.0f,		1.0f, 0.0f, 1.0f, c},
		{ 1.0f, 0.0f, 1.0f, 1.0f,    0.0f, 1.0f, 1.0f, 1.0f,    0.0f, 0.0f, 1.0f, 1.0f,		0.0f, 1.0f, 1.0f,		1.0f, 0.0f, 1.0f,		1.0f, 1.0f, 1.0f, c},
						   																			   
		// WEST            																			   
		{ 0.0f, 0.0f, 1.0f, 1.0f,    0.0f, 1.0f, 1.0f, 1.0f,    0.0f, 1.0f, 0.0f, 1.0f,		0.0f, 1.0f, 1.0f,		0.0f, 0.0f, 1.0f,		1.0f, 0.0f, 1.0f, c},
		{ 0.0f, 0.0f, 1.0f, 1.0f,    0.0f, 1.0f, 0.0f, 1.0f,    0.0f, 0.0f, 0.0f, 1.0f,		0.0f, 1.0f, 1.0f,		1.0f, 0.0f, 1.0f,		1.0f, 1.0f, 1.0f, c},
						   																			   
		// TOP             																			   
		{ 0.0f, 1.0f, 0.0f, 1.0f,    0.0f, 1.0f, 1.0f, 1.0f,    1.0f, 1.0f, 1.0f, 1.0f,		0.0f, 1.0f, 1.0f,		0.0f, 0.0f, 1.0f,		1.0f, 0.0f, 1.0f, c},
		{ 0.0f, 1.0f, 0.0f, 1.0f,    1.0f, 1.0f, 1.0f, 1.0f,    1.0f, 1.0f, 0.0f, 1.0f,		0.0f, 1.0f, 1.0f,		1.0f, 0.0f, 1.0f,		1.0f, 1.0f, 1.0f, c},
						   																			  
		// BOTTOM          																			  
		{ 1.0f, 0.0f, 1.0f, 1.0f,    0.0f, 0.0f, 1.0f, 1.0f,    0.0f, 0.0f, 0.0f, 1.0f,		0.0f, 1.0f, 1.0f,		0.0f, 0.0f, 1.0f,		1.0f, 0.0f, 1.0f, c},
		{ 1.0f, 0.0f, 1.0f, 1.0f,    0.0f, 0.0f, 0.0f, 1.0f,    1.0f, 0.0f, 0.0f, 1.0f,		0.0f, 1.0f, 1.0f,		1.0f, 0.0f, 1.0f,		1.0f, 1.0f, 1.0f, c},
	};

	for(int i = Plane*2;i<(Plane+1)*2;i++){
		for(int j = 0;j<3;j++){
			tris[i].p[j] = vec3d_Add(p,vec3d_Make(tris[i].p[j].x * d.x,tris[i].p[j].y * d.y,tris[i].p[j].z * d.z));
		}
		Vector_Push(&meshCube.tris,&tris[i]);
	}
}

void int_swap(int* a,int* b){
	int c = *a;
	*a = *b;
	*b = c;
}
void float_swap(float* a,float* b){
	float c = *a;
	*a = *b;
	*b = c;
}
float float_max(float a,float b){
	return a>b ? a : b;
}


void TexturedTriangle(	int x1, int y1, float u1, float v1, float w1,
						int x2, int y2, float u2, float v2, float w2,
						int x3, int y3, float u3, float v3, float w3,
                        unsigned char id,
						Sprite *tex){
	if (y2 < y1){
		int_swap(&y1,&y2);
		int_swap(&x1,&x2);
		float_swap(&u1,&u2);
		float_swap(&v1,&v2);
		float_swap(&w1,&w2);
	}
	if (y3 < y1){
		int_swap(&y1,&y3);
		int_swap(&x1,&x3);
		float_swap(&u1,&u3);
		float_swap(&v1,&v3);
		float_swap(&w1,&w3);
	}
	if (y3 < y2){
		int_swap(&y2,&y3);
		int_swap(&x2,&x3);
		float_swap(&u2,&u3);
		float_swap(&v2,&v3);
		float_swap(&w2,&w3);
	}

	int dy1 = y2 - y1;
	int dx1 = x2 - x1;
	float dv1 = v2 - v1;
	float du1 = u2 - u1;
	float dw1 = w2 - w1;
	int dy2 = y3 - y1;
	int dx2 = x3 - x1;
	float dv2 = v3 - v1;
	float du2 = u3 - u1;
	float dw2 = w3 - w1;
	float tex_u, tex_v, tex_w;
	float dax_step = 0, dbx_step = 0,
		du1_step = 0, dv1_step = 0,
		du2_step = 0, dv2_step = 0,
		dw1_step=0, dw2_step=0;
	if (dy1) dax_step = dx1 / (float)fabsf(dy1);
	if (dy2) dbx_step = dx2 / (float)fabsf(dy2);
	if (dy1) du1_step = du1 / (float)fabsf(dy1);
	if (dy1) dv1_step = dv1 / (float)fabsf(dy1);
	if (dy1) dw1_step = dw1 / (float)fabsf(dy1);
	if (dy2) du2_step = du2 / (float)fabsf(dy2);
	if (dy2) dv2_step = dv2 / (float)fabsf(dy2);
	if (dy2) dw2_step = dw2 / (float)fabsf(dy2);
	
	if (dy1){
		for (int i = y1; i <= y2; i++){
			int ax = x1 + (float)(i - y1) * dax_step;
			int bx = x1 + (float)(i - y1) * dbx_step;
			float tex_su = u1 + (float)(i - y1) * du1_step;
			float tex_sv = v1 + (float)(i - y1) * dv1_step;
			float tex_sw = w1 + (float)(i - y1) * dw1_step;
			float tex_eu = u1 + (float)(i - y1) * du2_step;
			float tex_ev = v1 + (float)(i - y1) * dv2_step;
			float tex_ew = w1 + (float)(i - y1) * dw2_step;
			
			if (ax > bx){
				int_swap(&ax,&bx);
				float_swap(&tex_su,&tex_eu);
				float_swap(&tex_sv,&tex_ev);
				float_swap(&tex_sw,&tex_ew);
			}
			
			tex_u = tex_su;
			tex_v = tex_sv;
			tex_w = tex_sw;
			float tstep = 1.0f / ((float)(bx - ax));
			float t = 0.0f;
			for (int j = ax; j < bx; j++){
				tex_u = (1.0f - t) * tex_su + t * tex_eu;
				tex_v = (1.0f - t) * tex_sv + t * tex_ev;
				tex_w = (1.0f - t) * tex_sw + t * tex_ew;
				if (tex_w > pDepthBuffer[i*GetWidth() + j]){
					//Draw(j,i,Sprite_Sample(tex,tex_u / tex_w, tex_v / tex_w));
					int subx = id % 16;
					int suby = id / 16;
					Draw(j,i,Sprite_SampleSub(tex,tex_u / tex_w, tex_v / tex_w,subx,suby,16,16));
					pDepthBuffer[i*GetWidth() + j] = tex_w;
				}
				t += tstep;
			}
		}
	}

	dy1 = y3 - y2;
	dx1 = x3 - x2;
	dv1 = v3 - v2;
	du1 = u3 - u2;
	dw1 = w3 - w2;
	if (dy1) dax_step = dx1 / (float)fabsf(dy1);
	if (dy2) dbx_step = dx2 / (float)fabsf(dy2);
	du1_step = 0, dv1_step = 0;
	if (dy1) du1_step = du1 / (float)fabsf(dy1);
	if (dy1) dv1_step = dv1 / (float)fabsf(dy1);
	if (dy1) dw1_step = dw1 / (float)fabsf(dy1);
	
	if(dy1){
		for(int i = y2; i <= y3; i++){
			int ax = x2 + (float)(i - y2) * dax_step;
			int bx = x1 + (float)(i - y1) * dbx_step;
			float tex_su = u2 + (float)(i - y2) * du1_step;
			float tex_sv = v2 + (float)(i - y2) * dv1_step;
			float tex_sw = w2 + (float)(i - y2) * dw1_step;
			float tex_eu = u1 + (float)(i - y1) * du2_step;
			float tex_ev = v1 + (float)(i - y1) * dv2_step;
			float tex_ew = w1 + (float)(i - y1) * dw2_step;
			
			if(ax > bx){
				int_swap(&ax,&bx);
				float_swap(&tex_su,&tex_eu);
				float_swap(&tex_sv,&tex_ev);
				float_swap(&tex_sw,&tex_ew);
			}

			tex_u = tex_su;
			tex_v = tex_sv;
			tex_w = tex_sw;
			float tstep = 1.0f / ((float)(bx - ax));
			float t = 0.0f;
			for(int j = ax; j < bx; j++){
				tex_u = (1.0f - t) * tex_su + t * tex_eu;
				tex_v = (1.0f - t) * tex_sv + t * tex_ev;
				tex_w = (1.0f - t) * tex_sw + t * tex_ew;
				if (tex_w > pDepthBuffer[i*GetWidth() + j])
				{
					int subx = id % 16;
					int suby = id / 16;
					Draw(j,i,Sprite_SampleSub(tex,tex_u / tex_w, tex_v / tex_w,subx,suby,16,16));
					pDepthBuffer[i * GetWidth() + j] = tex_w;
				}
				t += tstep;
			}
		}	
	}		
}

void Setup(AlxWindow* w){
	pDepthBuffer = malloc(sizeof(float) * GetWidth() * GetHeight());

	meshCube = mesh_New();

	MakeCube(vec3d_Make(0.0f,0.0f,0.0f),vec3d_Make(1.0f,1.0f,1.0f),WHITE);
	MakeCube(vec3d_Make(0.0f,0.0f,1.0f),vec3d_Make(1.0f,1.0f,1.0f),WHITE);
	MakeCube(vec3d_Make(1.0f,0.0f,0.0f),vec3d_Make(1.0f,1.0f,1.0f),WHITE);
	MakeCube(vec3d_Make(1.0f,0.0f,1.0f),vec3d_Make(1.0f,1.0f,1.0f),WHITE);

	vCamera = vec3d_Make(0.0f,0.0f,-4.0f);
	vLookDir = vec3d_New();
	fYaw = 0.0f;
	fTheta = 0.0f;
	
	sprTex1 = Sprite_Load("./Atlas.png");
	matProj = Matrix_MakeProjection(90.0f,(float)GetHeight() / (float)GetWidth(),0.1f,1000.0f);
}

void Update(AlxWindow* w){
	if (Stroke(ALX_KEY_UP).DOWN)
		vCamera.y += 8.0f * w->ElapsedTime;
	if (Stroke(ALX_KEY_DOWN).DOWN)
		vCamera.y -= 8.0f * w->ElapsedTime;
	
	vec3d vForward = vec3d_Mul(vLookDir,8.0f * w->ElapsedTime);
	if (Stroke(ALX_KEY_W).DOWN)
		vCamera = vec3d_Add(vCamera,vForward);
	if (Stroke(ALX_KEY_S).DOWN)
		vCamera = vec3d_Sub(vCamera,vForward);
	if (Stroke(ALX_KEY_A).DOWN)
		fYaw -= 2.0f * w->ElapsedTime;
	if (Stroke(ALX_KEY_D).DOWN)
		fYaw += 2.0f * w->ElapsedTime;
	
	fTheta += 1.0f * w->ElapsedTime;
	mat4x4 matRotZ = Matrix_MakeRotationZ(fTheta * 0.5f);
	mat4x4 matRotX = Matrix_MakeRotationX(fTheta);
	
	mat4x4 matTrans = Matrix_MakeTranslation(0.0f,0.0f,0.0f);
	
	mat4x4 matWorld = Matrix_MakeIdentity();
	matWorld = Matrix_MultiplyMatrix(matRotZ, matRotX);
	matWorld = Matrix_MultiplyMatrix(matWorld, matTrans);
	
	vec3d vUp = vec3d_Make(0.0f,1.0f,0.0f);
	vec3d vTarget = vec3d_Make(0.0f,0.0f,1.0f);
	
	mat4x4 matCameraRot = Matrix_MakeRotationY(fYaw);
	vLookDir = Matrix_MultiplyVector(matCameraRot,vTarget);
	
	mat4x4 matCamera = Matrix_PointAt(vCamera,vTarget,vUp);
	mat4x4 matView = Matrix_QuickInverse(matCamera);
	
	Vector vecTrianglesToRaster = Vector_New(sizeof(triangle));
	for(int i = 0;i<meshCube.tris.size;i++){
		triangle tri = *(triangle*)Vector_Get(&meshCube.tris,i);

		triangle triTransformed;
		triTransformed.p[0] = Matrix_MultiplyVector(matWorld, tri.p[0]);
		triTransformed.p[1] = Matrix_MultiplyVector(matWorld, tri.p[1]);
		triTransformed.p[2] = Matrix_MultiplyVector(matWorld, tri.p[2]);
		triTransformed.t[0] = tri.t[0];
		triTransformed.t[1] = tri.t[1];
		triTransformed.t[2] = tri.t[2];
		
		vec3d normal, line1, line2;
		line1 = vec3d_Sub(triTransformed.p[1], triTransformed.p[0]);
		line2 = vec3d_Sub(triTransformed.p[2], triTransformed.p[0]);
		normal = vec3d_CrossProduct(line1, line2);
		normal = vec3d_Normalise(normal);
		
		vec3d vCameraRay = vec3d_Sub(triTransformed.p[0], vCamera);
		if (vec3d_DotProduct(normal, vCameraRay) < 0.0f){
			vec3d light_direction = vec3d_Make(0.0f,1.0f,-1.0f);
			light_direction = vec3d_Normalise(light_direction);
			
			float dp = float_max(0.1f,vec3d_DotProduct(light_direction,normal));
			
			Pixel c = Pixel_toRGBA(dp,dp,dp,1.0f);
			triTransformed.c = c;
			triTransformed.id = tri.id;

			triangle triViewed;
			triViewed.p[0] = Matrix_MultiplyVector(matView,triTransformed.p[0]);
			triViewed.p[1] = Matrix_MultiplyVector(matView,triTransformed.p[1]);
			triViewed.p[2] = Matrix_MultiplyVector(matView,triTransformed.p[2]);
			triViewed.c = triTransformed.c;
			triViewed.id = triTransformed.id;
			triViewed.t[0] = triTransformed.t[0];
			triViewed.t[1] = triTransformed.t[1];
			triViewed.t[2] = triTransformed.t[2];
			
			
			int nClippedTriangles = 0;
			triangle clipped[2];
			nClippedTriangles = Triangle_ClipAgainstPlane(vec3d_Make(0.0f,0.0f,0.1f),vec3d_Make(0.0f,0.0f,1.0f),&triViewed,&clipped[0],&clipped[1]);
			
			for (int n = 0; n < nClippedTriangles; n++){
				triangle triProjected;
				triProjected.p[0] = Matrix_MultiplyVector(matProj,clipped[n].p[0]);
				triProjected.p[1] = Matrix_MultiplyVector(matProj,clipped[n].p[1]);
				triProjected.p[2] = Matrix_MultiplyVector(matProj,clipped[n].p[2]);
				triProjected.c = clipped[n].c;
				triProjected.id = clipped[n].id;
				triProjected.t[0] = clipped[n].t[0];
				triProjected.t[1] = clipped[n].t[1];
				triProjected.t[2] = clipped[n].t[2];

				triProjected.t[0].u = triProjected.t[0].u / triProjected.p[0].w;
				triProjected.t[1].u = triProjected.t[1].u / triProjected.p[1].w;
				triProjected.t[2].u = triProjected.t[2].u / triProjected.p[2].w;
				triProjected.t[0].v = triProjected.t[0].v / triProjected.p[0].w;
				triProjected.t[1].v = triProjected.t[1].v / triProjected.p[1].w;
				triProjected.t[2].v = triProjected.t[2].v / triProjected.p[2].w;
				triProjected.t[0].w = 1.0f / triProjected.p[0].w;
				triProjected.t[1].w = 1.0f / triProjected.p[1].w;
				triProjected.t[2].w = 1.0f / triProjected.p[2].w;
				
				triProjected.p[0] = vec3d_Div(triProjected.p[0],triProjected.p[0].w);
				triProjected.p[1] = vec3d_Div(triProjected.p[1],triProjected.p[1].w);
				triProjected.p[2] = vec3d_Div(triProjected.p[2],triProjected.p[2].w);
				
				triProjected.p[0].x *= -1.0f;
				triProjected.p[1].x *= -1.0f;
				triProjected.p[2].x *= -1.0f;
				triProjected.p[0].y *= -1.0f;
				triProjected.p[1].y *= -1.0f;
				triProjected.p[2].y *= -1.0f;
				
				vec3d vOffsetView = vec3d_Make(1.0f,1.0f,0.0f);
				triProjected.p[0] = vec3d_Add(triProjected.p[0],vOffsetView);
				triProjected.p[1] = vec3d_Add(triProjected.p[1],vOffsetView);
				triProjected.p[2] = vec3d_Add(triProjected.p[2],vOffsetView);
				triProjected.p[0].x *= 0.5f * (float)GetWidth();
				triProjected.p[0].y *= 0.5f * (float)GetHeight();
				triProjected.p[1].x *= 0.5f * (float)GetWidth();
				triProjected.p[1].y *= 0.5f * (float)GetHeight();
				triProjected.p[2].x *= 0.5f * (float)GetWidth();
				triProjected.p[2].y *= 0.5f * (float)GetHeight();

				Vector_Push(&vecTrianglesToRaster,&triProjected);
			}			
		}
	}
	
	Clear(BLACK);
	memset(pDepthBuffer,0,sizeof(float) * GetWidth() * GetHeight());
		
	for(int i = 0;i<vecTrianglesToRaster.size;i++){
		triangle triToRaster = *(triangle*)Vector_Get(&vecTrianglesToRaster,i);
		
		triangle clipped[2];
		Vector listTriangles = Vector_New(sizeof(triangle));
		
		Vector_Push(&listTriangles,&triToRaster);
		int nNewTriangles = 1;
		for (int p = 0; p < 4; p++){
			int nTrisToAdd = 0;
			while (nNewTriangles > 0){
				triangle test = *(triangle*)Vector_Get(&listTriangles,0);
				Vector_Remove(&listTriangles,0);
				nNewTriangles--;
				
				switch (p)
				{
				case 0:	nTrisToAdd = Triangle_ClipAgainstPlane((vec3d){ 0.0f, 0.0f, 0.0f }, 					(vec3d){ 0.0f, 1.0f, 0.0f },&test,&clipped[0],&clipped[1]); break;
				case 1:	nTrisToAdd = Triangle_ClipAgainstPlane((vec3d){ 0.0f, (float)GetHeight() - 1, 0.0f }, 	(vec3d){ 0.0f,-1.0f, 0.0f },&test,&clipped[0],&clipped[1]); break;
				case 2:	nTrisToAdd = Triangle_ClipAgainstPlane((vec3d){ 0.0f, 0.0f, 0.0f }, 					(vec3d){ 1.0f, 0.0f, 0.0f },&test,&clipped[0],&clipped[1]); break;
				case 3:	nTrisToAdd = Triangle_ClipAgainstPlane((vec3d){ (float)GetWidth() - 1, 0.0f, 0.0f }, 	(vec3d){-1.0f, 0.0f, 0.0f },&test,&clipped[0],&clipped[1]); break;
				}
				
				for (int w = 0; w < nTrisToAdd; w++)
					Vector_Push(&listTriangles,&clipped[w]);
			}
			nNewTriangles = listTriangles.size;
		}
		
		for (int i = 0;i<listTriangles.size;i++){
			triangle t = *(triangle*)Vector_Get(&listTriangles,i);
			
			TexturedTriangle(
				t.p[0].x, t.p[0].y, t.t[0].u, t.t[0].v, t.t[0].w,
				t.p[1].x, t.p[1].y, t.t[1].u, t.t[1].v, t.t[1].w,
				t.p[2].x, t.p[2].y, t.t[2].u, t.t[2].v, t.t[2].w, 
				t.id,&sprTex1
			);
			
			//RenderTriangle(((Vec2){t.p[0].x,t.p[0].y}),((Vec2){t.p[1].x,t.p[1].y}),((Vec2){t.p[2].x,t.p[2].y}),t.c);
			//RenderTriangleWire(((Vec2){t.p[0].x,t.p[0].y}),((Vec2){t.p[1].x,t.p[1].y}),((Vec2){t.p[2].x,t.p[2].y}),t.c,1.0f);
			//RenderTriangleWire(((Vec2){t.p[0].x,t.p[0].y}),((Vec2){t.p[1].x,t.p[1].y}),((Vec2){t.p[2].x,t.p[2].y}),WHITE,1.0f);
		}
		Vector_Free(&listTriangles);
	}
	Vector_Free(&vecTrianglesToRaster);
}

void Delete(AlxWindow* w){
	free(pDepthBuffer);

	mesh_Free(&meshCube);
	Sprite_Free(&sprTex1);
}

int main(){
    if(Create("3D Test Tex",2500,1200,1,1,Setup,Update,Delete))
        Start();
    return 0;
}