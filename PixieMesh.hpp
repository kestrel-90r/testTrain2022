﻿# pragma once

# include <omp.h>
# include <thread>

# include <Siv3D.hpp>   
# include <Siv3D/EngineLog.hpp>
# include "PixieCamera.hpp"

#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define JSON_NOEXCEPTION
#define STBI_MSC_SECURE_CRT
#include "3rd/tiny_gltf.h"

using namespace DirectX;
#define CPUSKINNING
#define BASIS 0			

using Word4 = Vector4D<uint16>;

constexpr float USE_USERCOLOR = 0;	
constexpr float USE_OFFSET_METARIAL = -1;
constexpr float USE_TEXTURE = -2 ;

enum MODELTYPE { MODELNOA, MODELANI, MODELVRM };

enum Use {
	USE = -2, NOTUSE = -1,
	UNKNOWN = 0,
	USE_MORPH = 1, NOTUSE_MORPH,
	USE_VFP, NOTUSE_VFP,
	USE_INDIVIDUAL, NOTUSE_INDIVIDUAL,
	SHOW_COLLIDER, HIDE_COLLIDER, NOTUSE_COLLIDER
};


struct Channel
{
    uint8 typeDelta=0;
	uint8 numMorph=0;
    int16 idxNode=0;
    int32 idxSampler=0;
    Array<std::pair<int32, Array<float> >> deltaKeyframes;
	Channel() { numMorph = 0; typeDelta = 0; idxSampler = -1; idxNode = -1; }
};

struct Sampler
{
    Array<std::pair<int32, float>> interpolateKeyframes;
    float  minTime=0.0;
    float  maxTime=0.0;
    Sampler() { minTime = 0.0; maxTime = 0.0; }
};

struct Frame
{
    Array<MeshData>		MeshDatas;			 
    Array<DynamicMesh>	Meshes;				 
    Array<uint8>		useTex;				 
    Array<Mat4x4>		morphMatBuffers;	 
    Float3				obSize{1,1,1};		 
    Float3				obCenter{0,0,0};	 
};

struct NodeParam	
{
	Mat4x4 matLocal;								 
	Mat4x4 matWorld;								 
	Float3 posePos{0,0,0};							 
	Quaternion poseRot{ Quaternion::Identity() };	 
	Float3 poseSca{1,1,1};							 
	Mat4x4 matModify;								 
	Float3 modPos{0,0,0};							 
	Quaternion modRot{ Quaternion::Identity() };	 
	Float3 modSca{0,0,0};							 
	bool update=false;								 
};

struct PrecAnime      
{
    Array<ColorF>			meshColors;  
    Array<Texture>			meshTexs;    
    Array<Frame>			Frames;      

	Array<Array<Sampler>>   Samplers;    
    Array<Array<Channel>>   Channels;    
};

struct MorphMesh
{
    Array<int32>			Targets;        
	Array<Array<Vertex3D>>	ShapeBuffers;   		
    Array<Array<Vertex3D>>	BasisBuffers;   	
    uint32                   TexCoordCount; 
    Array<Float2>            TexCoord;	    
};

struct MorphTarget
{
	float Weight=0;					 
};

struct AnimeModel
{
    Array<PrecAnime>        precAnimes;
    MorphMesh               morphMesh;
};

struct NoAModel
{
    Array<String>           meshName;		 
    Array<ColorF>           meshColors;		 
    Array<Texture>          meshTexs;		 
    Array<MeshData>			MeshDatas;		 
    Array<DynamicMesh>		Meshes;			 
    Array<int32>            useTex;			 

    Array<Mat4x4>           morphMatBuffers; 
    MorphMesh               morphMesh;
};

struct VRMModel
{
    Array<String>           meshName;    
    Array<ColorF>           meshColors;  
    Array<Texture>          meshTexs;    
    Array<MeshData>			MeshDatas;	 
    Array<DynamicMesh>		Meshes;		 
    Array<int32>            useTex;      

    Array<Array<Mat4x4>>    Joints;
    Array<Mat4x4>           morphMatBuffers; 
    MorphMesh               morphMesh;
};

Quaternion QPYR(Float3 pyr = { 0,0,0 })
{
	return Quaternion::RollPitchYaw(ToRadians(pyr.x), ToRadians(pyr.y), ToRadians(pyr.z));
}
Quaternion QPYR(float p = 0, float y = 0, float r = 0)
{
	return Quaternion::RollPitchYaw(ToRadians(p), ToRadians(y), ToRadians(r));
}

#define DISPLACEMENTFUNC void (*displacementFunc)( Array<Vertex3D> &vertices, Array<TriangleIndex32> &indices )
class PixieMesh
{
private:
	NoAModel    noaModel;
	AnimeModel  aniModel;
	VRMModel    vrmModel;
	bool		register1st = false;

	Array<NodeParam> nodeParams;
	DISPLACEMENTFUNC = nullptr;
public:
	PixieCamera camera;  

	tinygltf::Model gltfModel;
	Array<MorphTarget>	morphTarget;	


	Use			showCollider = NOTUSE_COLLIDER;
	Use			effectDisplace = NOTUSE_VFP;

    String		gltfFilename = U"";

    Float3		Pos{0,0,0};						
    Float3		Sca{1,1,1};

	Float3		eRot{0,0,0};					
	Quaternion	qRot{ Quaternion::Identity() };	  
	Quaternion	qSpin{ Quaternion::Identity() };

    Float3		rPos{0,0,0};					

    int32		currentFrame=0;

	Array<OrientedBox>  obCollider;

	Mat4x4		Mat = Mat4x4::Identity();
	Mat4x4		matVP = Mat4x4::Identity();

    int32		animeID = 0;
    int32		morphID = 0;

	virtual ~PixieMesh() = default;
    PixieMesh() = default;
	explicit PixieMesh( String gltfFilename,   
			   Float3 pos = Float3{0,0,0},	   
               Float3 scale  = Float3{1,1,1},  
               Float3 erot = Float3{0,0,0},    
			   Float3 rotate = Float3{ 0,0,0 },
			   Float3 relpos = Float3{0,0,0},  
			   Float3 qspin = Float3{ 0,0,0 }, 
			   int32 frame=0,                  
               int32 anime=0,                  
               int32 morph=0)                  
    {
        gltfFilename = gltfFilename;
        Pos = pos;
        rPos = relpos;
        Sca = scale;
        eRot = erot;
		qRot = Quaternion::RollPitchYaw(ToRadians(rotate.x), ToRadians(rotate.y), ToRadians(rotate.z));

		if (qspin != Float3{ 0,0,0 })
			qSpin = Quaternion::RollPitchYaw(ToRadians(qspin.x), ToRadians(qspin.y), ToRadians(qspin.z));

		currentFrame = frame;
        animeID = anime;
        morphID = morph;

		camera.setEyePos(Pos + rPos);
		camera.setFocusPos(Pos + rPos + Float3{0,0,1});
	}

	PixieMesh& applyMove(Float3 amount)
	{
		Pos += amount;
		return *this;
	}
	PixieMesh& applyMoveRelative(Float3 amount)
	{
		rPos += amount;
		return *this;
	}
	PixieMesh& applyScale(Float3 amount)
	{
		Sca += amount;
		return *this;
	}
	PixieMesh& applyRotateEuler(Float3 amount)
	{
		eRot += amount;
		return *this;
	}
	PixieMesh& applyRotateQ(Quaternion amount)
	{
		qRot *= amount;
		return *this;
	}
	PixieMesh& applyMat4x4(Mat4x4 mat)
	{
		Float3 sca, tra;
		Quaternion rot;
		mat.decompose(sca, rot, tra);
		Pos += tra;
		qRot *= rot;
		Sca *= sca;
		return *this;
	}

	PixieMesh& setMove(Float3 amount)
	{
		Pos = amount;
		return *this;
	}
	PixieMesh& setMoveRelative(Float3 amount)
	{
		rPos = amount;
		return *this;
	}
	PixieMesh& setScale(Float3 amount)
	{
		Sca = amount;
		return *this;
	}
	PixieMesh& setRotateEuler(Float3 amount)
	{
		eRot = amount;
		return *this;
	}
	PixieMesh& setRotateQ(Quaternion amount)
	{
		qRot = amount;
		return *this;
	}
	PixieMesh& setSpinQ(Quaternion amount)
	{
		qSpin = amount;
		return *this;
	}
	PixieMesh& setShowCollider(Use use)
	{
		showCollider = use;
		return *this;
	}
	Use getShowCollisionBox()
	{
		return showCollider;
	}

	PixieMesh& setDisplacementFunc(DISPLACEMENTFUNC = nullptr)
	{
		this->displacementFunc = displacementFunc;
		return *this;
	}

	void initModel( String filename, MODELTYPE modeltype, const Size& sceneSize,
		Use useindividual = NOTUSE_INDIVIDUAL, Use usemorph = NOTUSE_MORPH,
		uint32 cycleframe = 60, int32 animeid = -1)
	{
		gltfFilename = filename;
		initModel(modeltype, sceneSize, useindividual, usemorph, cycleframe, animeid);
	}

    void initModel( MODELTYPE modeltype, const Size& sceneSize,
		Use useindividual = NOTUSE_INDIVIDUAL, Use usemorph = NOTUSE_MORPH, 
		uint32 cycleframe = 60, int32 animeid=-1)
    {
        std::string err, warn;
		bool result;
		tinygltf::TinyGLTF loader;

		result = loader.LoadBinaryFromFile(&gltfModel, &err, &warn, gltfFilename.narrow());

		if (!result) result = loader.LoadASCIIFromFile(&gltfModel, &err, &warn, gltfFilename.narrow());

		if (usemorph == NOTUSE_MORPH) morphTarget.clear();
		else
		{
			uint32 nmorph = 0;
			for (auto& mesh : gltfModel.meshes) nmorph += mesh.weights.size();
			morphTarget.resize(nmorph);
		}

		if (result && modeltype == MODELNOA)	  gltfSetupNOA(useindividual);
		else if (result && modeltype == MODELANI) gltfSetupANI( cycleframe, animeid);
		else if (result && modeltype == MODELVRM) gltfSetupVRM();

		camera = PixieCamera(sceneSize, 45_deg, Pos+rPos, Pos + rPos+Float3{ 0,0,1 }, 0.05);
    }

    void setStartFrame( uint32 anime_no, int32 offsetframe )
    {
        PrecAnime &pa = aniModel.precAnimes[anime_no];
        currentFrame = offsetframe % pa.Frames.size() ;
    }

    tinygltf::Buffer* getBuffer(tinygltf::Model& gltfmodel, tinygltf::Primitive& pr, size_t* offset, int* stride, int* componenttype)
    {
        if (pr.indices == -1) return nullptr;
        tinygltf::Accessor& ai = gltfmodel.accessors[pr.indices];
        tinygltf::BufferView& bi = gltfmodel.bufferViews[ai.bufferView];
        tinygltf::Buffer& buf = gltfmodel.buffers[bi.buffer];
        *offset = bi.byteOffset + ai.byteOffset;
        *stride = ai.ByteStride(bi);
        *componenttype = ai.componentType;
        return &buf;
    }

    tinygltf::Buffer* getBuffer(tinygltf::Model& gltfmodel, tinygltf::Primitive& pr, const std::string motion, size_t* offset, int* stride, int* componenttype)
    {
        if (pr.attributes.size() == 0) return nullptr;
        tinygltf::Accessor& ap = gltfmodel.accessors[pr.attributes[motion]];
        tinygltf::BufferView& bp = gltfmodel.bufferViews[ap.bufferView];
        tinygltf::Buffer& buf = gltfmodel.buffers[bp.buffer];
        *offset = bp.byteOffset + ap.byteOffset;
        *stride = ap.ByteStride(bp);
        *componenttype = ap.componentType;
        return &buf;
    }

    tinygltf::Buffer* getBuffer(tinygltf::Model& gltfmodel, tinygltf::Primitive& pr, const std::string motion, size_t* offset, int* stride, int* componenttype, int VARTYPE )
    {
        if (pr.attributes.size() == 0) return nullptr;
        tinygltf::Accessor& ap = gltfmodel.accessors[ pr.attributes[motion] ];
        tinygltf::BufferView& bp = gltfmodel.bufferViews[ ap.bufferView ];

		tinygltf::Buffer& buf = gltfmodel.buffers[bp.buffer];
        *offset = bp.byteOffset + ap.byteOffset;
		*componenttype = ap.componentType;
		*stride = ap.ByteStride(bp) ? (ap.ByteStride(bp) / tinygltf::GetComponentSizeInBytes(ap.componentType)) : tinygltf::GetNumComponentsInType(VARTYPE);
        return &buf;
    }

    static tinygltf::Buffer* getBuffer(tinygltf::Model& model, tinygltf::Primitive& pr, const int32 &morphtarget, const char* motion, size_t* offset, int* stride, int* componenttype)
    {
        if (pr.targets.size() == 0) return nullptr;
        tinygltf::Accessor& ap = model.accessors[pr.targets[morphtarget].at(motion)];
        tinygltf::BufferView& bp = model.bufferViews[ap.bufferView];
        tinygltf::Buffer& buf = model.buffers[bp.buffer];
        *offset = bp.byteOffset + ap.byteOffset;
        *stride = ap.ByteStride(bp);
        *componenttype = ap.componentType;
        return &buf;
    }

	void gltfSetupPosture(tinygltf::Model& gltfmodel, int32 nodeidx, Array<NodeParam>& nodeParams, Use useindividual = NOTUSE_INDIVIDUAL)
    {
		NodeParam &np = nodeParams[nodeidx];

		auto& node = gltfmodel.nodes[nodeidx];
		auto& r = node.rotation;
        auto& t = node.translation;
        auto& s = node.scale;

		np.matModify = np.matLocal = Mat4x4::Identity();
        Quaternion rr = Quaternion::Identity();
		Float3     tt = Float3{ 0,0,0 };
		Float3     ss = Float3{ 1,1,1 };

		if (useindividual != USE_INDIVIDUAL)
		{
			if (node.rotation.size()) rr = Quaternion{ r[0], r[1], r[2], r[3] };
			if (node.translation.size()) tt = Float3{ t[0], t[1], t[2] };
			if (node.scale.size()) ss = Float3{ s[0], s[1], s[2] };

			np.matLocal = Mat4x4::Identity().Scale(ss) *
						  Mat4x4(Quaternion(rr)) *
						  Mat4x4::Identity().Translate(tt);
		}

		np.matWorld = Mat4x4::Identity();	 
		np.posePos = Float3{0,0,0};			 
		np.poseRot = Quaternion::Identity(); 
		np.poseSca = Float3{1,1,1};			 

		np.update = false;					 

		for (uint32 cc = 0; cc < node.children.size(); cc++)
            gltfSetupPosture(gltfmodel, node.children[cc], nodeParams);
	}

    void gltfCalcSkeleton(tinygltf::Model& gltfmodel, const Mat4x4& matparent,
		                  int32 nodeidx,Array<NodeParam>& nodeParams )
    {
		NodeParam &np = nodeParams[nodeidx];

		auto& node = gltfmodel.nodes[nodeidx];
		Mat4x4 matlocal = np.matLocal ;		  

		Quaternion rr = np.poseRot;			  
		Float3 tt = np.posePos;
		Float3 ss = np.poseSca;
		Mat4x4 matpose = Mat4x4::Identity().Scale(ss)*
			             Mat4x4(rr)*
			             Mat4x4::Identity().Translate(tt);

		if (np.update)
		{
			matlocal = matlocal * np.matModify;
		}

		if( matpose.isIdentity() )
			matpose = matlocal;

		Mat4x4 mat = matpose * matlocal.inverse();
        Mat4x4 matworld = mat * matlocal * matparent;
		np.matWorld = matworld;

		for (uint32 cc = 0; cc < node.children.size(); cc++)
            gltfCalcSkeleton(gltfmodel, matworld, node.children[cc], nodeParams );
	}

    PixieMesh& gltfSetupNOA( Use useindividual = NOTUSE_INDIVIDUAL)
    {
		nodeParams.resize( gltfModel.nodes.size() );

		for (uint32 nn = 0; nn < gltfModel.nodes.size(); nn++)
			gltfSetupPosture( gltfModel, nn, nodeParams, useindividual);

		if (useindividual != USE_INDIVIDUAL)
		{
			for (uint32 nn = 0; nn < gltfModel.nodes.size(); nn++)
			{
				auto& mn = gltfModel.nodes[nn];
				if (mn.mesh >= 0)
					gltfCalcSkeleton(gltfModel, Mat4x4::Identity(), nn, nodeParams);

				for (uint32 cc = 0; cc < mn.children.size(); cc++)
					gltfCalcSkeleton(gltfModel, Mat4x4::Identity(), mn.children[cc], nodeParams);
			}
		}

		Array<Array<Mat4x4>> Joints;
		if (gltfModel.skins.size() > 0)
		{
			Joints.resize(gltfModel.skins.size());
			for (uint32 nn = 0; nn < gltfModel.nodes.size(); nn++)
			{
				auto& mn = gltfModel.nodes[nn];
				if (mn.skin >= 0)         
				{
					Joints[ mn.skin ].resize( gltfModel.skins[mn.skin].joints.size() );

					auto& msns = gltfModel.skins[mn.skin];
					auto& ibma = gltfModel.accessors[msns.inverseBindMatrices];	
					auto& ibmbv = gltfModel.bufferViews[ibma.bufferView];		
					auto  ibmd = gltfModel.buffers[ibmbv.buffer].data.data() + ibma.byteOffset + ibmbv.byteOffset;

					for (uint32 ii = 0; ii < msns.joints.size(); ii++)
					{
						Mat4x4 ibm = *(Mat4x4*)&ibmd[ii * sizeof(Mat4x4)];
						Mat4x4& matworld = nodeParams[msns.joints[ii]].matWorld;
						Joints[mn.skin][ii] = ibm * matworld;

					}
				}
			}
		}


		for (uint32 nn = 0; nn < gltfModel.nodes.size(); nn++)
			gltfSetupMorph(gltfModel.nodes[nn], noaModel.morphMesh);

		for (uint32 nn = 0; nn < gltfModel.nodes.size(); nn++)
			gltfSetupNOA(gltfModel.nodes[nn], nn, Joints);

		obCollider.clear();
		if (useindividual != USE_INDIVIDUAL)
		{
			Float3 vmin = { FLT_MAX,FLT_MAX,FLT_MAX };
			Float3 vmax = { FLT_MIN,FLT_MIN,FLT_MIN };

			for (uint32 i = 0; i < noaModel.MeshDatas.size(); i++)
			{
				for (uint32 ii = 0; ii < noaModel.MeshDatas[i].vertices.size(); ii++)
				{
					Vertex3D& mv = noaModel.MeshDatas[i].vertices[ii];

					if (vmin.x > mv.pos.x) vmin.x = mv.pos.x;
					if (vmin.y > mv.pos.y) vmin.y = mv.pos.y;
					if (vmin.z > mv.pos.z) vmin.z = mv.pos.z;
					if (vmax.x < mv.pos.x) vmax.x = mv.pos.x;
					if (vmax.y < mv.pos.y) vmax.y = mv.pos.y;
					if (vmax.z < mv.pos.z) vmax.z = mv.pos.z;
				}
			}

			obCollider.emplace_back();
			obCollider[0].setOrientation(QPYR(eRot) * qRot);
			obCollider[0].setPos( vmin + (vmax - vmin) / 2);
			obCollider[0].setSize(vmax - vmin);

		}

		else
		{
			for (uint32 i = 0; i < noaModel.MeshDatas.size(); i++)
			{
				obCollider.emplace_back( noaModel.Meshes[i].boundingBox() );
			}
		}

		nodeParams.clear();
		return *this;
	}

    void gltfSetupNOA( tinygltf::Node& node, uint32 nodeidx, Array<Array<Mat4x4>> &Joints )
    {
		static uint32 morphidx = 0;         

		if (node.mesh < 0) return;			
		else if (node.mesh == 0)
		{
			morphidx = 0;
		}

		auto& weights = gltfModel.meshes[node.mesh].weights ;
		uint32 prsize = (uint32)gltfModel.meshes[node.mesh].primitives.size();
        for (uint32 pp = 0; pp < prsize; pp++)
        {
            auto& pr = gltfModel.meshes[node.mesh].primitives[pp];

			size_t opos=0, otex=0, onor=0, ojoints=0, oweights=0, oidx=0;
            int32 type_p=0, type_t=0, type_n=0, type_j=0, type_w=0, type_i=0;
			int32 stride_p=0, stride_t=0, stride_n=0, stride_j=0, stride_w=0, stride_i=0;

			auto& bpos = *getBuffer(gltfModel, pr, "POSITION", &opos, &stride_p, &type_p);        
			auto& btex = *getBuffer(gltfModel, pr, "TEXCOORD_0", &otex, &stride_t, &type_t);	  
			auto& bnor = *getBuffer(gltfModel, pr, "NORMAL", &onor, &stride_n, &type_n);       
			auto& bjoint = *getBuffer(gltfModel, pr, "JOINTS_0", &ojoints, &stride_j, &type_j, TINYGLTF_TYPE_VEC4);   
			auto& bweight = *getBuffer(gltfModel, pr, "WEIGHTS_0", &oweights, &stride_w, &type_w);
			auto& bidx = *getBuffer(gltfModel, pr, &oidx, &stride_i, &type_i);                    


            Array<Vertex3D> vertices;
			auto& numvertex = gltfModel.accessors[pr.attributes["POSITION"]].count; 
			for (uint32 vv = 0; vv < numvertex; vv++)
            {
				Vertex3D mv;
				float* pos = (float*)&bpos.data.at(vv * stride_p + opos);
				float* nor = (float*)&bnor.data.at(vv * stride_n + onor);
				float* tex = (float*)&btex.data.at(vv * stride_t + otex);

				mv.pos = Float3{ pos[0], pos[1], pos[2] };  
				mv.normal = Float3{ nor[0], nor[1], nor[2] };
				mv.tex = Float2{ tex[0], tex[1] };


				if ( pr.targets.size() )
				{
					for (uint32 tt = 0; tt < pr.targets.size(); tt++)
					{
						if (weights[tt] == 0.0) continue;

						size_t opos=0, onor=0;
						auto& mtpos = *getBuffer(gltfModel, pr, tt, "POSITION", &opos, &stride_p, &type_p);
						auto& mtnor = *getBuffer(gltfModel, pr, tt, "NORMAL", &onor, &stride_n, &type_n);
						float* spos = (float*)&mtpos.data.at(vv * stride_p + opos);
						float* snor = (float*)&mtnor.data.at(vv * stride_n + onor);
						Float3 shapepos = Float3(spos[0], spos[1], spos[2]);
						Float3 shapenor = Float3(snor[0], snor[1], snor[2]);
						mv.pos += shapepos * weights[tt];
						mv.normal += shapenor * weights[tt];
					}
				}

				Mat4x4 matlocal = nodeParams[nodeidx].matLocal;
                SIMD_Float4 vec4pos = DirectX::XMVector4Transform(SIMD_Float4(mv.pos, 1.0f), matlocal);
                Mat4x4 matnor = matlocal.inverse().transposed();
                mv.normal = SIMD_Float4{ DirectX::XMVector4Transform(SIMD_Float4(mv.normal, 1.0f), matlocal) }.xyz();

				Mat4x4 matskin = Mat4x4::Identity();
				if( node.skin >= 0 )
                {

#define USE_GLTF
#ifdef USE_GLTF
					uint8* jb = (uint8*)&bjoint.data.at(vv * stride_j + ojoints); 
					uint16* jw = (uint16*)&bjoint.data.at(vv * stride_j + ojoints);
					Word4 j4 = (type_j == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT) ?
											Word4(jw[0], jw[1], jw[2], jw[3]) :
											Word4(jb[0], jb[1], jb[2], jb[3]);
#endif

#ifdef USE_VRM
					float* jd = (float*)&bjoint.data.at(vv * stride_j + ojoints);
					Word4 j4 = Word4(jd[0],jd[1],jd[2],jd[3]) ;
#endif

					float* wf = (float*)&bweight.data.at(vv * stride_w + oweights);
					Float4 w4 = Float4(wf[0], wf[1], wf[2], wf[3]);

                    matskin = w4.x * Joints[node.skin][j4.x] +
						   	  w4.y * Joints[node.skin][j4.y] +
							  w4.z * Joints[node.skin][j4.z] +
							  w4.w * Joints[node.skin][j4.w];

					if (pr.targets.size() > 0)
						noaModel.morphMatBuffers.emplace_back(matskin);
					vec4pos = DirectX::XMVector4Transform(SIMD_Float4(mv.pos, 1.0f), matskin);
					if (!(std::abs(matskin.determinant()) < 10e-10))
					matnor = matskin.inverse().transposed();
                }

				mv.pos = vec4pos.xyz() /vec4pos.getW();
				mv.normal = SIMD_Float4{ DirectX::XMVector4Transform(SIMD_Float4(mv.normal, 1.0f), matskin) }.xyz();
                vertices.emplace_back(mv);
            }

            MeshData md;
            if (pr.indices > 0)
            {
                auto& mapi = gltfModel.accessors[pr.indices];
                Array<TriangleIndex32> indices;
                for (int32 i = 0; i < mapi.count; i+=3)
                {
					TriangleIndex32 idx = TriangleIndex32::Zero();
					if (mapi.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT)
					{
						uint16* ibuf = (uint16 *)&bidx.data.at(i * 2 + oidx); 
						idx.i0 = ibuf[0]; idx.i1 = ibuf[1]; idx.i2 = ibuf[2];
					}
					else if (mapi.componentType == TINYGLTF_COMPONENT_TYPE_UNSIGNED_INT)
					{
						uint32* ibuf = (uint32 *)&bidx.data.at(i * 4 + oidx); 
						idx.i0 = ibuf[0]; idx.i1 = ibuf[1]; idx.i2 = ibuf[2];
					}
                    indices.emplace_back(idx);
                }

                md = MeshData(vertices, indices);
                vertices.clear();
                indices.clear();
            }

            int32 usetex = 0;   
			Texture tex ;
            ColorF col = ColorF(1);

            if (pr.material >= 0)
            {
                auto& nt = gltfModel.materials[pr.material].additionalValues["normalTexture"];  
                int32 idx = -1;
                auto& mmv = gltfModel.materials[pr.material].values;
                auto& bcf = mmv["baseColorFactor"];												

				if (mmv.count("baseColorTexture"))
                    idx = gltfModel.textures[(int32)mmv["baseColorTexture"].json_double_value["index"]].source;

                if (bcf.number_array.size()) col = ColorF(  bcf.number_array[0],
                                                            bcf.number_array[1],
                                                            bcf.number_array[2],
                                                            bcf.number_array[3]);
                else                         col = ColorF(1);

                if (idx >= 0 && gltfModel.images.size())
                {
                    tex = Texture();
                    if (gltfModel.images[idx].bufferView >= 0)	 
                    {
                        auto& bgfx = gltfModel.bufferViews[gltfModel.images[idx].bufferView];
                        auto bimg = &gltfModel.buffers[bgfx.buffer].data.at(bgfx.byteOffset);

						tex = Texture( MemoryReader { bimg,bgfx.byteLength}, TextureDesc::MippedSRGB);
						usetex = 1;
					}

					else
                    {
                    }

                }
            }

			noaModel.meshTexs.emplace_back(tex);								 
			noaModel.meshColors.emplace_back(col);                               
            noaModel.meshName.emplace_back(Unicode::FromUTF8(gltfModel.meshes[node.mesh].name));   
            noaModel.MeshDatas.emplace_back(md);                                 
			noaModel.Meshes.emplace_back( DynamicMesh{ md } );                   
            noaModel.useTex.emplace_back( usetex );                              
        }
	}

    void gltfSetupMorph( tinygltf::Node& node, MorphMesh& morph )
    {
		if (node.mesh < 0) return;		
		uint32 prsize = (uint32)gltfModel.meshes[node.mesh].primitives.size();
		for (uint32 pp = 0; pp < prsize; pp++)
        {
            auto& pr = gltfModel.meshes[node.mesh].primitives[pp];
			morph.Targets.emplace_back((int32)pr.targets.size());	

			if (pr.targets.size() > 0)	
            {
				size_t opos = 0, otex = 0, onor = 0;
				int32 stride_p = 0, stride_t = 0, stride_n = 0;
				int32 type_p = 0, type_t=0, type_n = 0;

                auto bpos = getBuffer(gltfModel, pr, "POSITION", &opos, &stride_p, &type_p);
                auto bnor = getBuffer(gltfModel, pr, "NORMAL", &onor, &stride_n, &type_n);
				auto btex = getBuffer(gltfModel, pr, "TEXCOORD_0", &otex, &stride_t, &type_t);	  


                Array<Vertex3D> basisvertices;
                auto& numvertex = gltfModel.accessors[pr.attributes["POSITION"]].count; 
                for (int32 vv = 0; vv < numvertex; vv++)
                {
                    Vertex3D mv;
					float* pos = (float*)&bpos->data.at(vv * stride_p + opos);
					float* nor = (float*)&bnor->data.at(vv * stride_n + onor);
					float* tex = (float*)&btex->data.at(vv * stride_t + otex);
					mv.pos = Float3(pos[0], pos[1], pos[2]);
                    mv.normal = Float3(nor[0], nor[1], nor[2]);
					mv.tex = Float2(tex[0], tex[1]);
                    basisvertices.emplace_back(mv);

				}

                morph.BasisBuffers.emplace_back(basisvertices);
                basisvertices.clear();

				for (int32 tt = 0; tt < pr.targets.size(); tt++)
                {
					size_t opos = 0, otex = 0, onor = 0;
					int32 stride=0;
                    int32 type_p=0, type_n=0;

                    auto mtpos = getBuffer(gltfModel, pr, tt, "POSITION", &opos, &stride, &type_p);
                    auto mtnor = getBuffer(gltfModel, pr, tt, "NORMAL", &onor, &stride, &type_n);

					Array<Vertex3D> shapevertices;
                    for (int32 vv = 0; vv < numvertex; vv++)
                    {
						Vertex3D mv;
                        auto pos = (float*)&mtpos->data.at(vv * 12 + opos);
                        auto nor = (float*)&mtnor->data.at(vv * 12 + onor);
						auto tex = (float*)&btex->data.at(vv * 8 + otex);
						mv.pos = Float3(pos[0], pos[1], pos[2]);
                        mv.normal = Float3(nor[0], nor[1], nor[2]);
						mv.tex = Float2(tex[0], tex[1]);
						shapevertices.emplace_back(mv);
                    }
                    morph.ShapeBuffers.emplace_back(shapevertices);
                    shapevertices.clear();
                }

            }
        }
    }

	void gltfSetupVRM()
	{
		nodeParams.resize( gltfModel.nodes.size() );

		for (int32 nn = 0; nn < gltfModel.nodes.size(); nn++)
			gltfSetupPosture( gltfModel, nn, nodeParams );

		for (int32 nn = 0; nn < gltfModel.nodes.size(); nn++)
		{
			auto& mn = gltfModel.nodes[nn];

			Mat4x4 imat = Mat4x4::Identity();
			if (mn.mesh >= 0)
				gltfCalcSkeleton(gltfModel, imat, nn, nodeParams );

			for (int32 cc = 0; cc < mn.children.size(); cc++)
				gltfCalcSkeleton(gltfModel, imat, mn.children[cc], nodeParams  );
		}

		if (gltfModel.skins.size() > 0)
		{
			vrmModel.Joints.resize(gltfModel.skins.size());
			for (uint32 nn = 0; nn < gltfModel.nodes.size(); nn++)
			{
				auto& mn = gltfModel.nodes[nn];
				if (mn.skin >= 0)         
				{
					vrmModel.Joints[mn.skin].resize(gltfModel.skins[mn.skin].joints.size());

					const auto& msns = gltfModel.skins[mn.skin];
					const auto& ibma = gltfModel.accessors[msns.inverseBindMatrices];
					const auto& ibmbv = gltfModel.bufferViews[ibma.bufferView];
					auto  ibmd = &gltfModel.buffers[ibmbv.buffer].data.data()[ ibma.byteOffset + ibmbv.byteOffset];

					for (uint32 ii = 0; ii < msns.joints.size(); ii++)
					{
						Mat4x4 ibm = *(Mat4x4*)&ibmd[ii * sizeof(Mat4x4)];
						Mat4x4& matworld = nodeParams[msns.joints[ii]].matWorld;
						Mat4x4 matjoint = ibm * matworld;

						vrmModel.Joints[mn.skin][ii] = matjoint;
					}
				}
			}
		}

		for (uint32 nn = 0; nn < gltfModel.nodes.size(); nn++)
			gltfSetupMorph(gltfModel.nodes[nn], vrmModel.morphMesh);

		for (uint32 nn = 0; nn < gltfModel.nodes.size(); nn++)
			gltfSetupVRM(gltfModel.nodes[nn], nn);


		register1st = true;	
	}

	PixieMesh &drawVRM(uint32 istart = -1, uint32 icount = -1)
	{
		if (Pos.hasNaN() || qRot.hasNaN() || qRot.hasInf()) return *this;

		for (int32 nn = 0; nn < gltfModel.nodes.size(); nn++)
        {
            auto& mn = gltfModel.nodes[ nn ];

			Mat4x4 imat = Mat4x4::Identity();
            if (mn.mesh >= 0)
                gltfCalcSkeleton(gltfModel, imat, nn, nodeParams );

            for (int32 cc = 0; cc < mn.children.size(); cc++)
                gltfCalcSkeleton(gltfModel, imat, mn.children[cc], nodeParams );
        }

		for (int32 nn = 0; nn < gltfModel.nodes.size(); nn++)
        {
            auto& node = gltfModel.nodes[ nn ];

            if (node.skin < 0) continue;         

            const auto& msns = gltfModel.skins[node.skin];
            const auto& ibma = gltfModel.accessors[msns.inverseBindMatrices];
            const auto& ibmbv = gltfModel.bufferViews[ibma.bufferView];
            auto  ibmd = gltfModel.buffers[ibmbv.buffer].data.data() + ibma.byteOffset + ibmbv.byteOffset;

			for (int32 ii = 0; ii < vrmModel.Joints[node.skin].size(); ii++)
            {
                Mat4x4 ibm = *(Mat4x4*)&ibmd[ii * sizeof(Mat4x4)];
				Mat4x4 &matworld = nodeParams[ msns.joints[ii] ].matWorld;
				Mat4x4 matjoint = ibm * matworld ;

				vrmModel.Joints[node.skin][ii] = matjoint;
			}
        }

		for (uint32 nn = 0; nn < gltfModel.nodes.size(); nn++)
		{
            auto& mn = gltfModel.nodes[nn];
			if (mn.mesh >= 0)
            {
                gltfSetupVRM( mn, nn );
				gltfDrawVRM( istart, icount );
            }
        }
		return *this;
    }

	void gltfSetupVRM(tinygltf::Node& node, uint32 nodeidx )
	{
		static uint32 morphidx = 0;         
		static uint32 meshidx = 0;			
		if (node.mesh < 0) return;
		else if (node.mesh == 0)
		{
			morphidx = 0;
			meshidx = 0;
		}

		Array<Vertex3D> morphmv;
		const Array<Array<Vertex3D>>& shapes = vrmModel.morphMesh.ShapeBuffers;

		int32 prsize = gltfModel.meshes[node.mesh].primitives.size();
		for (int32 pp = 0; pp < prsize; pp++)
		{
			auto& pr = gltfModel.meshes[node.mesh].primitives[pp];
			auto& map = gltfModel.accessors[pr.attributes["POSITION"]];

			size_t opos=0, otex=0, onormal=0, ojoints=0, oweights=0, oidx=0;
			int32 stride_p = 0, stride_n = 0, stride_t = 0, stride_j = 0, stride_w = 0,stride = 0;
			int32 type_p=0, type_t=0, type_n=0, type_j=0, type_w=0, type_i=0;

			auto bpos = getBuffer(gltfModel, pr, "POSITION", &opos, &stride_p, &type_p);        
			auto btex = getBuffer(gltfModel, pr, "TEXCOORD_0", &otex, &stride_t, &type_t);	  
			auto bnormal = getBuffer(gltfModel, pr, "NORMAL", &onormal, &stride_n, &type_n);       
			auto bjoint = getBuffer(gltfModel, pr, "JOINTS_0", &ojoints, &stride_j, &type_j);   
			auto bweight = getBuffer(gltfModel, pr, "WEIGHTS_0", &oweights, &stride_w, &type_w);
			auto bidx = getBuffer(gltfModel, pr, &oidx, &stride, &type_i);                    

			const uint32 NMORPH = pr.targets.size();
			if (NMORPH > 0 )  
			{
				morphmv = vrmModel.morphMesh.BasisBuffers[morphidx];

				for (uint32 vv = 0; vv < morphmv.size(); vv++)	 
				{
					for (uint32 mm = 0; mm < NMORPH; mm++)			 
					{
						float& weight = morphTarget[mm].Weight;
						if (weight == 0.0) continue;
						morphmv[vv].pos += shapes[pp * NMORPH + mm][vv].pos * weight;
						morphmv[vv].normal += shapes[pp * NMORPH + mm][vv].normal * weight;
					}
				}
                morphidx++;
			}

			Array<Vertex3D> vertices;
			for (int32 vv = 0; vv < map.count; vv++)
			{
				float* vertex = nullptr, * texcoord = nullptr, * normal = nullptr;
				vertex = (float*)&bpos->data.at(vv * stride_p + opos);
				texcoord = (float*)&btex->data.at(vv * stride_t + otex);
				normal = (float*)&bnormal->data.at(vv * stride_n + onormal);


				Vertex3D mv;
				if ( pr.targets.size() > 0) mv = morphmv[vv];									 
				else						mv.pos = Float3{ vertex[0], vertex[1], vertex[2] };  

				mv.tex = Float2{ texcoord[0], texcoord[1] };
				mv.normal = Float3{ normal[0], normal[1], normal[2] };

				Mat4x4 matlocal = nodeParams[nodeidx].matLocal;
                SIMD_Float4 vec4pos = DirectX::XMVector4Transform(SIMD_Float4(mv.pos, 1.0f), matlocal);
                Mat4x4 matnor = matlocal.inverse().transposed();
                mv.normal = SIMD_Float4{ DirectX::XMVector4Transform(SIMD_Float4(mv.normal, 1.0f), matlocal) }.xyz();

				Mat4x4 matskin = Mat4x4::Identity();
				if( node.skin >= 0 )
				{
					uint8* jb = (uint8*)&bjoint->data.at(vv * 4 + ojoints); 
					uint16* jw = (uint16*)&bjoint->data.at(vv * 8 + ojoints);
					Word4 j4 = (type_j == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT) ? Word4(jw[0], jw[1], jw[2], jw[3]) :  Word4(jb[0], jb[1], jb[2], jb[3]) ;

					float* wf = (float*)&bweight->data.at(vv * 16 + oweights);
					Float4 w4 = Float4(wf[0], wf[1], wf[2], wf[3]);

					matskin = w4.x * vrmModel.Joints[node.skin][j4.x] +
						      w4.y * vrmModel.Joints[node.skin][j4.y] +
						      w4.z * vrmModel.Joints[node.skin][j4.z] +
						      w4.w * vrmModel.Joints[node.skin][j4.w];

					if (pr.targets.size() > 0)
					{
						if (register1st == false) vrmModel.morphMatBuffers.emplace_back(matskin);
						else					  vrmModel.morphMatBuffers[pp] = matskin;
					}
					vec4pos = DirectX::XMVector4Transform(SIMD_Float4(mv.pos, 1.0f), matskin);
					matnor = matskin.inverse().transposed();
				}

				mv.pos = Float3(vec4pos.getX(), vec4pos.getY(), vec4pos.getZ()) / vec4pos.getW();
				mv.normal = SIMD_Float4{ DirectX::XMVector4Transform(SIMD_Float4(mv.normal, 1.0f), matskin) }.xyz();
				vertices.emplace_back(mv);
            }

			if (register1st == false)
			{
				MeshData md;
				auto& mapi = gltfModel.accessors[pr.indices];
				Array<TriangleIndex32> indices;
				for (int32 i = 0; i < mapi.count; i += 3)
				{
					TriangleIndex32 idx = TriangleIndex32::Zero();
					if (mapi.componentType == 5123)
					{
						uint16* ibuf = (uint16*)&bidx->data.at(i * 2 + oidx); 
						idx.i0 = ibuf[0]; idx.i1 = ibuf[1]; idx.i2 = ibuf[2];
					}
					else if (mapi.componentType == 5125)
					{
						uint32* ibuf = (uint32*)&bidx->data.at(i * 4 + oidx); 
						idx.i0 = ibuf[0]; idx.i1 = ibuf[1]; idx.i2 = ibuf[2];
					}
					indices.emplace_back(idx);
				}
				md = MeshData(vertices, indices);
				indices.clear();
				vertices.clear();
				vrmModel.MeshDatas.emplace_back(md);				 
				vrmModel.Meshes.emplace_back( DynamicMesh{ md });	 
			}

			else
			{
				vrmModel.Meshes[meshidx++].fill(vertices);
				vertices.clear();
			}

			if ( register1st == false )
			{
				int32 usetex = 0;   
				Texture tex;
				ColorF col = ColorF(1, 1, 1, 1);
				if (pr.material >= 0)
				{
					auto& nt = gltfModel.materials[pr.material].additionalValues["normalTexture"];  
					auto& mmv = gltfModel.materials[pr.material].values;
					auto& bcf = mmv["baseColorFactor"];												
					int32 idx = -1;
					if (mmv.count("baseColorTexture"))
						idx = gltfModel.textures[(int32)mmv["baseColorTexture"].json_double_value["index"]].source;

					if (bcf.number_array.size()) col = ColorF(bcf.number_array[0],
															  bcf.number_array[1],
															  bcf.number_array[2],
															  bcf.number_array[3]);

					if (idx >= 0 && gltfModel.images.size())
					{
						tex = Texture();
						if (gltfModel.images[idx].bufferView >= 0)	 
						{
							auto& bgfx = gltfModel.bufferViews[gltfModel.images[idx].bufferView];
							auto bimg = &gltfModel.buffers[bgfx.buffer].data.at(bgfx.byteOffset);
							tex = Texture(MemoryReader{ bimg,bgfx.byteLength }, TextureDesc::MippedSRGB);
							usetex = 1;
						}

						else
						{
						}

					}
				}

				vrmModel.meshTexs.emplace_back(tex);	 
				vrmModel.meshColors.emplace_back(col);	 
				vrmModel.meshName.emplace_back(Unicode::FromUTF8(gltfModel.meshes[node.mesh].name));   
				vrmModel.useTex.emplace_back(usetex);      
			}
		}
	}

	void gltfDrawVRM( uint32 istart = -1, uint32 icount = -1 )
	{
        uint32 tid = 0;
        __m128 rot = XMQuaternionRotationRollPitchYaw(ToRadians(eRot.x), ToRadians(eRot.y), ToRadians(eRot.z));

		Mat4x4 mrot;
		if (!qSpin.isIdentity()) mrot = Mat4x4(Quaternion(rot) * qRot * qSpin);
		else mrot = Mat4x4(Quaternion(rot) * qRot);

		Float3 tra = Pos + rPos;

		Mat = Mat4x4::Identity().Scale(Float3{ -Sca.x,Sca.y,Sca.z }) * mrot *
		      Mat4x4::Identity().Translate(tra);

		for (uint32 i = 0; i < vrmModel.Meshes.size(); i++)
        {
			if (istart == -1)						
			{
				if (vrmModel.useTex[i])  vrmModel.Meshes[i].draw(Mat, vrmModel.meshTexs[i], vrmModel.meshColors[i]);	
				else					 vrmModel.Meshes[i].draw(Mat, vrmModel.meshColors[i]);							
			}
			else									
			{
				if (vrmModel.useTex[i])  vrmModel.Meshes[i].drawSubset(istart, icount, Mat, vrmModel.meshTexs[tid++]);	
				else					 vrmModel.Meshes[i].drawSubset(istart, icount, Mat, vrmModel.meshColors[i]);	
			}
        }
    }


    PixieMesh &setCamera(const PixieCamera& cam)
    {
        camera = cam;
        matVP = cam.getViewProj();
        return *this;
    }


	static void drawGridBox(Mat4x4 mat = Mat4x4::Identity(), int32 size = 10000, String text = U"")
	{
		float s = size;
		for (int32 i = 0; i < s; i+=size/10)
		{
			Float3 xy[2] = { mat.transformPoint({-s,i,0}) , mat.transformPoint({+s,i,0}) };
			Float3 xz[2] = { mat.transformPoint({-s,0,i}) , mat.transformPoint({+s,0,i}) };
			Float3 yx[2] = { mat.transformPoint({i,-s,0}) , mat.transformPoint({i,+s,0}) };
			Float3 yz[2] = { mat.transformPoint({0,-s,i}) , mat.transformPoint({0,+s,i}) };
			Float3 zx[2] = { mat.transformPoint({i,0,-s}) , mat.transformPoint({i,0,+s}) };
			Float3 zy[2] = { mat.transformPoint({0,i,-s}) , mat.transformPoint({0,i,+s}) };
			Line3D(xy[0], xy[1]).draw(ColorF{ 0.8,0.2,0.2 });
			Line3D(xz[0], xz[1]).draw(ColorF{ 0.8,0.2,0.2 });
			Line3D(yx[0], yx[1]).draw(ColorF{ 0.2,0.8,0.2 });
			Line3D(yz[0], yz[1]).draw(ColorF{ 0.2,0.8,0.2 });
			Line3D(zx[0], zx[1]).draw(ColorF{ 0.2,0.2,0.8 });
			Line3D(zy[0], zy[1]).draw(ColorF{ 0.2,0.2,0.8 });
		}
	}

	PixieMesh& setMat(Mat4x4 &mat)
	{
		Float3 pos, sca;
		Quaternion qrot;
		mat.decompose(sca,qrot,pos);
		setRotateQ(qrot);
		setMove(pos);
		setScale(sca);
		return *this;
	}

	const Mat4x4 &getMat()
	{
		return Mat;
	}

	PixieMesh& drawMesh(const Mat4x4 mat, int32 meshcode = NOTUSE,
		ColorF usrColor = ColorF(0, 0, 0, USE_TEXTURE), int32 istart = NOTUSE, int32 icount = NOTUSE)
	{
		if (Pos.hasNaN() || qRot.hasNaN() || qRot.hasInf()) return *this;

		NoAModel& noa = noaModel;
		uint32 primitiveidx = 0;
		uint32 tid = 0;

		Array<Vertex3D> morphmv;
		const Array<Array<Vertex3D>>& shapes = noa.morphMesh.ShapeBuffers;

		const uint32 NMORPH = morphTarget.size();

		Mat = mat;

		uint32 i = meshcode;
		if (meshcode == NOTUSE || meshcode >= noa.Meshes.size())
		{
			for (i = 0; i < noa.Meshes.size(); i++)
			{
#if 0
				if (noa.morphMesh.Targets[i])
				{
					morphmv = noa.morphMesh.BasisBuffers[primitiveidx]; 

					for (uint32 vv = 0; vv < morphmv.size(); vv++)	
					{
						for (uint32 mm = 0; mm < morphTarget.size(); mm++)   
						{
							float& weight = morphTarget[mm].Weight;
							if (weight == 0.0) continue;
							morphmv[vv].pos += shapes[primitiveidx * NMORPH + mm][vv].pos * weight;
							morphmv[vv].normal += shapes[primitiveidx * NMORPH + mm][vv].normal * weight;
						}

						Mat4x4 matskin = Mat4x4::Identity();
						Mat4x4 matnor = Mat4x4::Identity();
						SIMD_Float4 vec4pos = DirectX::XMVector4Transform(SIMD_Float4(morphmv[vv].pos, 1.0f), matskin);
						if (noa.morphMatBuffers.size())
						{
							matskin = noa.morphMatBuffers[vv];
							matnor = matskin.inverse().transposed();
						}

						morphmv[vv].pos = vec4pos.xyz() / vec4pos.getW();
						morphmv[vv].normal = SIMD_Float4{ DirectX::XMVector4Transform(SIMD_Float4(morphmv[vv].normal, 1.0f), matnor) }.xyz();
						morphmv[vv].tex = noa.MeshDatas[i].vertices[vv].tex;
					}

					noa.Meshes[i].fill(morphmv);				
					primitiveidx++;

				}
#endif
#if 0
				if (displacementFunc != nullptr)
				{
					Array<Vertex3D>	vertices = noa.MeshDatas[i].vertices;	
					Array<TriangleIndex32> indices = noa.MeshDatas[i].indices;
					(*displacementFunc)(vertices, indices);
					noa.Meshes[i].fill(vertices);							
					noa.Meshes[i].fill(indices);
				}
#endif
				if (istart == NOTUSE)	
				{
					if (usrColor.a >= USE_USERCOLOR)
						noa.Meshes[i].draw(Mat, ColorF(usrColor.rgb(), 1));

					else
					{
						if (usrColor.a == USE_OFFSET_METARIAL)	
							noa.Meshes[i].draw(Mat, noa.meshColors[i] + ColorF(usrColor.rgb(), 1));
						else
							noa.Meshes[i].draw(Mat, noa.meshTexs[i], noa.meshColors[i]);
					}
				}
				else					
				{
					if (usrColor.a >= USE_USERCOLOR)
						noa.Meshes[i].drawSubset(istart, icount, Mat, ColorF(usrColor.rgb(), 1));
					else
					{
						if (usrColor.a == USE_OFFSET_METARIAL)	
							noa.Meshes[i].drawSubset(istart, icount, Mat, noa.meshColors[i] + ColorF(usrColor.rgb(), 1));
						else
							noa.Meshes[i].drawSubset(istart, icount, Mat, noa.meshTexs[i], noa.meshColors[i]);
					}
				}


			}
		}

		else
		{
#if 0
			if (noa.morphMesh.Targets[i])
			{
				morphmv = noa.morphMesh.BasisBuffers[primitiveidx]; 

				for (uint32 vv = 0; vv < morphmv.size(); vv++)	
				{
					for (uint32 mm = 0; mm < morphTarget.size(); mm++)   
					{
						float& weight = morphTarget[mm].Weight;
						if (weight == 0.0) continue;
						morphmv[vv].pos += shapes[primitiveidx * NMORPH + mm][vv].pos * weight;
						morphmv[vv].normal += shapes[primitiveidx * NMORPH + mm][vv].normal * weight;
					}

					Mat4x4 matskin = Mat4x4::Identity();
					Mat4x4 matnor = Mat4x4::Identity();
					SIMD_Float4 vec4pos = DirectX::XMVector4Transform(SIMD_Float4(morphmv[vv].pos, 1.0f), matskin);
					if (noa.morphMatBuffers.size())
					{
						matskin = noa.morphMatBuffers[vv];
						matnor = matskin.inverse().transposed();
					}

					morphmv[vv].pos = vec4pos.xyz() / vec4pos.getW();
					morphmv[vv].normal = SIMD_Float4{ DirectX::XMVector4Transform(SIMD_Float4(morphmv[vv].normal, 1.0f), matnor) }.xyz();
					morphmv[vv].tex = noa.MeshDatas[i].vertices[vv].tex;
				}

				noa.Meshes[i].fill(morphmv);				
				primitiveidx++;

			}
#endif
#if 0
			if (displacementFunc != nullptr)
			{
				Array<Vertex3D>	vertices = noa.MeshDatas[i].vertices;	
				Array<TriangleIndex32> indices = noa.MeshDatas[i].indices;
				(*displacementFunc)(vertices, indices);
				noa.Meshes[i].fill(vertices);							
				noa.Meshes[i].fill(indices);
			}
#endif
			if (istart == NOTUSE)	
			{
				if (usrColor.a >= USE_USERCOLOR)
					noa.Meshes[i].draw(Mat, noa.meshTexs[i], noa.meshColors[i]);

				else
				{
					if (usrColor.a == USE_OFFSET_METARIAL)	
						noa.Meshes[i].draw(Mat, noa.meshColors[i] + ColorF(usrColor.rgb(), 1));
					else
						noa.Meshes[i].draw(Mat, noa.meshTexs[i], noa.meshColors[i]);
				}
			}
			else				
			{
				if (usrColor.a >= USE_USERCOLOR)
					noa.Meshes[i].drawSubset(istart, icount, Mat, ColorF(usrColor.rgb(), 1));
				else				
				{
					if (usrColor.a == USE_OFFSET_METARIAL)	
						noa.Meshes[i].drawSubset(istart, icount, Mat, noa.meshColors[i] + ColorF(usrColor.rgb(), 1));
					else
						noa.Meshes[i].drawSubset(istart, icount, Mat, noa.meshTexs[i], noa.meshColors[i]);
				}
			}
		}


		i = 0;
		if ( meshcode != NOTUSE ) i = meshcode;

		OrientedBox box = obCollider[i];
		if (showCollider != NOTUSE_COLLIDER)
			box = Geometry3D::TransformBoundingOrientedBox(obCollider[i], Mat);

		if (showCollider == SHOW_COLLIDER)
			box.drawFrame(Palette::Blue);


		return *this;
	}

	PixieMesh& drawMesh(int32 meshcode = NOTUSE, ColorF usrColor = ColorF(0,0,0,USE_TEXTURE), int32 istart = NOTUSE, int32 icount = NOTUSE)
	{
		if (Pos.hasNaN() || qRot.hasNaN() || qRot.hasInf()) return *this;

		Rect rectdraw = Rect{ 0,0,camera.getSceneSize() };

		NoAModel& noa = noaModel;
		uint32 primitiveidx = 0;
		uint32 tid = 0;

		if (!qSpin.isIdentity()) qRot *= qSpin;
		Quaternion rot = qRot ;
		if (eRot != Float3{0,0,0}) rot *= QPYR(eRot) ;

		Mat = Mat4x4::Identity().Scale(Sca) * Mat4x4(rot) *
		      Mat4x4::Identity().Translate( Pos + rPos );

		Array<Vertex3D> morphmv;
		const Array<Array<Vertex3D>>& shapes = noa.morphMesh.ShapeBuffers;

		const uint32 NMORPH = morphTarget.size();

		if (meshcode == NOTUSE || meshcode >= noa.Meshes.size() )
		{
			for (uint32 i = 0; i < noa.Meshes.size(); i++)
			{
#if 0
				if (noa.morphMesh.Targets[i])
				{
					morphmv = noa.morphMesh.BasisBuffers[primitiveidx]; 

					for (uint32 vv = 0; vv < morphmv.size(); vv++)	
					{
						for (uint32 mm = 0; mm < morphTarget.size(); mm++)   
						{
							float& weight = morphTarget[mm].Weight;
							if (weight == 0.0) continue;
							morphmv[vv].pos += shapes[primitiveidx * NMORPH + mm][vv].pos * weight;
							morphmv[vv].normal += shapes[primitiveidx * NMORPH + mm][vv].normal * weight;
						}

						Mat4x4 matskin = Mat4x4::Identity();
						Mat4x4 matnor = Mat4x4::Identity();
						SIMD_Float4 vec4pos = DirectX::XMVector4Transform(SIMD_Float4(morphmv[vv].pos, 1.0f), matskin);
						if (noa.morphMatBuffers.size())
						{
							matskin = noa.morphMatBuffers[vv];
							matnor = matskin.inverse().transposed();
						}

						morphmv[vv].pos = vec4pos.xyz() / vec4pos.getW();
						morphmv[vv].normal = SIMD_Float4{ DirectX::XMVector4Transform(SIMD_Float4(morphmv[vv].normal, 1.0f), matnor) }.xyz();
						morphmv[vv].tex = noa.MeshDatas[i].vertices[vv].tex;
					}

					noa.Meshes[i].fill(morphmv);				
					primitiveidx++;

				}
#endif
#if 0
				if (displacementFunc != nullptr)
				{
					Array<Vertex3D>	vertices = noa.MeshDatas[i].vertices;	
					Array<TriangleIndex32> indices = noa.MeshDatas[i].indices;
					(*displacementFunc)(vertices, indices);
					noa.Meshes[i].fill(vertices);							
					noa.Meshes[i].fill(indices);
				}
#endif
				if (istart == NOTUSE)	
				{
					if (usrColor.a >= USE_USERCOLOR)		
						noa.Meshes[i].draw(Mat, ColorF(usrColor.rgb(), 1));

					else
					{
						if (usrColor.a == USE_OFFSET_METARIAL)	
							noa.Meshes[i].draw(Mat, noa.meshColors[i] + ColorF(usrColor.rgb(), 1));
						else
							noa.Meshes[i].draw(Mat, noa.meshTexs[i], noa.meshColors[i]);
					}
				}
				else				
				{
					if (usrColor.a >= USE_USERCOLOR)		
						noa.Meshes[i].drawSubset(istart, icount, Mat, ColorF(usrColor.rgb(), 1));
					else				
					{
						if (usrColor.a == USE_OFFSET_METARIAL)	
							noa.Meshes[i].drawSubset(istart, icount, Mat, noa.meshColors[i] + ColorF(usrColor.rgb(), 1));
						else
							noa.Meshes[i].drawSubset(istart, icount, Mat, noa.meshTexs[i], noa.meshColors[i]);
					}
				}
			}
		}

		else
		{
#if 0
			if (noa.morphMesh.Targets[i])
			{
				morphmv = noa.morphMesh.BasisBuffers[primitiveidx]; 

				for (uint32 vv = 0; vv < morphmv.size(); vv++)	
				{
					for (uint32 mm = 0; mm < morphTarget.size(); mm++)   
					{
						float& weight = morphTarget[mm].Weight;
						if (weight == 0.0) continue;
						morphmv[vv].pos += shapes[primitiveidx * NMORPH + mm][vv].pos * weight;
						morphmv[vv].normal += shapes[primitiveidx * NMORPH + mm][vv].normal * weight;
					}

					Mat4x4 matskin = Mat4x4::Identity();
					Mat4x4 matnor = Mat4x4::Identity();
					SIMD_Float4 vec4pos = DirectX::XMVector4Transform(SIMD_Float4(morphmv[vv].pos, 1.0f), matskin);
					if (noa.morphMatBuffers.size())
					{
						matskin = noa.morphMatBuffers[vv];
						matnor = matskin.inverse().transposed();
					}

					morphmv[vv].pos = vec4pos.xyz() / vec4pos.getW();
					morphmv[vv].normal = SIMD_Float4{ DirectX::XMVector4Transform(SIMD_Float4(morphmv[vv].normal, 1.0f), matnor) }.xyz();
					morphmv[vv].tex = noa.MeshDatas[i].vertices[vv].tex;
				}

				noa.Meshes[i].fill(morphmv);				
				primitiveidx++;

			}
#endif
#if 0
			if (displacementFunc != nullptr)
			{
				Array<Vertex3D>	vertices = noa.MeshDatas[i].vertices;	
				Array<TriangleIndex32> indices = noa.MeshDatas[i].indices;
				(*displacementFunc)(vertices, indices);
				noa.Meshes[i].fill(vertices);							
				noa.Meshes[i].fill(indices);
			}
#endif
			int32 i = meshcode;

			if (istart == NOTUSE)	
			{
				if ( usrColor.a >= USE_USERCOLOR )
					noa.Meshes[i].draw(Mat, ColorF(usrColor.rgb(), 1));

				else
				{
					if (usrColor.a == USE_OFFSET_METARIAL)	
						noa.Meshes[i].draw(Mat, noa.meshColors[i] + ColorF(usrColor.rgb(), 1));
					else
						noa.Meshes[i].draw(Mat, noa.meshTexs[i], noa.meshColors[i]);
				}
			}
			else				
			{

				if (usrColor.a >= USE_USERCOLOR)
					noa.Meshes[i].drawSubset(istart, icount, Mat, ColorF(usrColor.rgb(), 1));
				else
				{
					if (usrColor.a == USE_OFFSET_METARIAL)	
						noa.Meshes[i].drawSubset(istart, icount, Mat, noa.meshColors[i] + ColorF(usrColor.rgb(), 1));
					else
						noa.Meshes[i].drawSubset(istart, icount, Mat, noa.meshTexs[i], noa.meshColors[i]);
				}
			}
		}

		int i = 0;
		if ( meshcode != NOTUSE ) i = meshcode;

		OrientedBox box = obCollider[i];
		if (showCollider != NOTUSE_COLLIDER)
			box = Geometry3D::TransformBoundingOrientedBox(obCollider[i], Mat);

		if (showCollider == SHOW_COLLIDER)
			box.drawFrame(Palette::Blue);

		return *this;
	}

	PixieMesh& gltfSetupANI( uint32 cycleframe = 60, int32 animeid=-1 )
    {
		aniModel.precAnimes.resize(gltfModel.animations.size());
		aniModel.morphMesh.TexCoordCount = (unsigned)-1;

        if (animeid == -1)                      
        {
            for (int32 aid = 0; aid < gltfModel.animations.size(); aid++)
				gltfOmpSetupANI( aid, cycleframe);
        }
        else
            gltfOmpSetupANI( animeid, cycleframe); 

		return *this;
	}

	const OrientedBox &getCollider( int32 meshcode )
	{
		if ( meshcode == NOTUSE ) meshcode = 0;
		if ( meshcode > obCollider.size()) assert(1);
		return obCollider[meshcode];
	}

	Box getBoundBox(int32 meshcode )
	{
		if (meshcode == NOTUSE) meshcode = 0;
		if (meshcode > noaModel.Meshes.size()) assert(1);
		return noaModel.Meshes[meshcode].boundingBox();
	}

	PixieMesh& setCollider( int32 meshcode, const OrientedBox &collider )
	{
		if (meshcode == NOTUSE) meshcode = 0;
		if (meshcode > obCollider.size()) assert(1);
		obCollider[ meshcode ] = collider;
		return *this;
	}

	void gltfOmpSetupANI(int32 animeid = 0, int32 cycleframe = 60)
	{
        if (gltfModel.animations.size() == 0) return ;

		auto& gm = gltfModel;
		auto& man = gltfModel.animations[ animeid ];
        auto& mas = man.samplers;
        auto& mac = man.channels;

        auto& macc = gltfModel.accessors;
        auto& mbv = gltfModel.bufferViews;
        auto& mb = gltfModel.buffers;

        auto begintime = macc[mas[0].input].minValues[0];     
        auto endtime = macc[mas[0].input].maxValues[0];       
        auto frametime = (endtime - begintime) / cycleframe;  
        auto currenttime = begintime;

		aniModel.precAnimes[ animeid ].Frames.resize(cycleframe);

		Array<double> frametimes;
		for (int32 cf = 0; cf < cycleframe; cf++)
		{
			frametimes.emplace_back(currenttime += frametime);
		}

		omp_set_num_threads(32);
		int32 tmax = omp_get_max_threads();
		aniModel.precAnimes[ animeid ].Samplers.resize(tmax);
		aniModel.precAnimes[ animeid ].Channels.resize(tmax);

		size_t bufsize = 0;
		for (int32 th = 0; th < tmax; th++)
		{
			auto& as = aniModel.precAnimes[ animeid ].Samplers[th];
			auto& ac = aniModel.precAnimes[ animeid ].Channels[th];
			as.resize(mas.size()); bufsize += mas.size()*sizeof(as);
			ac.resize(mac.size()); bufsize += mac.size()*sizeof(ac);

			for (int32 ss = 0; ss < mas.size(); ss++)
			{
				auto& mbsi = gltfModel.accessors[mas[ss].input];   
				as[ss].interpolateKeyframes.resize(mbsi.count); bufsize += mbsi.count*sizeof(as[ss].interpolateKeyframes);
			}

			for (int32 cc = 0; cc < ac.size(); cc++)
			{
				auto& maso = gltfModel.accessors[mas[cc].output];
				ac[cc].deltaKeyframes.resize(maso.count);bufsize += maso.count *sizeof(ac[cc].deltaKeyframes);

				auto& mid = gltfModel.nodes[mac[cc].target_node].mesh;
				ac[cc].numMorph = 0;
				if (mid != -1) ac[cc].numMorph = (uint8)gltfModel.meshes[ mid ].weights.size();

				for (int32 ff = 0; ff < maso.count; ff++)
				{
					if (ac[cc].numMorph)
					{
						ac[cc].deltaKeyframes[ff].second.resize(ac[cc].numMorph);
						bufsize += ac[cc].numMorph * sizeof(ac[cc].deltaKeyframes[ff].second);
					}
					else
					{
						ac[cc].deltaKeyframes[ff].second.resize(4);
						bufsize += 4;
					}
				}
			}
			for (int32 ss = 0; ss < mas.size(); ss++)
			{
				auto& mssi = gltfModel.accessors[mas[ss].input];   
				as[ss].interpolateKeyframes.resize(mssi.count);
				bufsize += mssi.count * sizeof(as[ss].interpolateKeyframes);
			}

			if (th == 1 && aniModel.morphMesh.TexCoordCount == (unsigned)-1)
				aniModel.morphMesh.TexCoordCount = aniModel.morphMesh.TexCoord.size();
		}

#pragma omp parallel for
		for ( int32 cf = 0; cf < cycleframe; cf++)
		{
			int32 th = omp_get_thread_num();

			auto& frametime = frametimes[cf];

			auto& man = gm.animations[ animeid ];
			auto& mas = man.samplers;
			auto& mac = man.channels;

			auto& macc = gm.accessors;
			auto& mbv = gm.bufferViews;
			auto& mb = gm.buffers;

			auto& as = aniModel.precAnimes[ animeid ].Samplers[th];
			auto& ac = aniModel.precAnimes[ animeid ].Channels[th];

			Array<NodeParam> nodeAniParams( gm.nodes.size() );

			for (int32 nn = 0; nn < gm.nodes.size(); nn++)
				gltfSetupPosture(gm, nn, nodeAniParams );

			for (int32 ss = 0; ss < mas.size(); ss++)
			{
				auto& msi = gm.accessors[mas[ss].input];   

				as[ss].minTime = 0;
				as[ss].maxTime = 1;
				if (msi.minValues.size() > 0 && msi.maxValues.size() > 0)
				{
					as[ss].minTime = float(msi.minValues[0]);
					as[ss].maxTime = float(msi.maxValues[0]);
				}

				for (int32 kk = 0; kk < msi.count; kk++)
				{
					auto& bsi = mas[ss].input;
					const auto& offset = mbv[macc[bsi].bufferView].byteOffset + macc[bsi].byteOffset;
					const auto& stride = msi.ByteStride(mbv[msi.bufferView]);
					void* adr = &mb[mbv[macc[bsi].bufferView].buffer].data.at(offset + kk * stride);

					auto& ctype = macc[bsi].componentType;
					float value =   (ctype == 5126) ? *(float*)adr :
									(ctype == 5123) ? *(uint16*)adr :
									(ctype == 5121) ? *(uint8_t*)adr :
									(ctype == 5122) ? *(int16*)adr :
									(ctype == 5120) ? *(int8_t*)adr : 0.0;

					as[ss].interpolateKeyframes[kk] = std::make_pair(kk, value);
				}
			}

			for (int32 cc = 0; cc < ac.size(); cc++)
			{
				auto& maso = gm.accessors[mas[cc].output];
				auto& macso = macc[mas[mac[cc].sampler].output];
				const auto& stride = maso.ByteStride(mbv[maso.bufferView]);
				const auto& offset = mbv[macso.bufferView].byteOffset + macso.byteOffset;

				ac[cc].idxNode = (uint16)mac[cc].target_node;
				ac[cc].idxSampler = mac[cc].sampler;

				if ( gm.nodes[mac[cc].target_node].mesh != -1)
					ac[cc].numMorph = (uint8)gm.meshes[ gm.nodes[ mac[cc].target_node ].mesh ].weights.size() ;

				if (mac[cc].target_path == "weights")
				{
					ac[cc].typeDelta = 5;    

					for (int32 ff = 0; ff<maso.count/ac[cc].numMorph; ff++)      
					{
						float* weight = (float*)&mb[mbv[macso.bufferView].buffer].data.at(offset + ff * stride * ac[cc].numMorph);
						ac[cc].deltaKeyframes[ff].first = ff;

						for (int32 mm = 0; mm<ac[cc].numMorph; mm++)     
							ac[cc].deltaKeyframes[ff].second[mm] = weight[mm] ;
					}
				}

				else if (mac[cc].target_path == "translation")
				{
					ac[cc].typeDelta = 1;   

					for (int32 ff = 0; ff < maso.count; ff++)
					{
						float* tra = (float*)&mb[mbv[macso.bufferView].buffer].data.at(offset + ff * stride);
						ac[cc].deltaKeyframes[ff].first = ff;
						ac[cc].deltaKeyframes[ff].second[0] = tra[0];
						ac[cc].deltaKeyframes[ff].second[1] = tra[1];
						ac[cc].deltaKeyframes[ff].second[2] = tra[2];
					}
				}

				else if (mac[cc].target_path == "rotation")
				{
					ac[cc].typeDelta = 3;

					for (int32 ff = 0; ff < maso.count; ff++)
					{
						float* rot = (float*)&mb[mbv[macso.bufferView].buffer].data.at(offset + ff * stride);
						auto qt = Quaternion(rot[0], rot[1], rot[2], rot[3]).normalize();

						ac[cc].deltaKeyframes[ff].first = ff;
						ac[cc].deltaKeyframes[ff].second[0] = qt.getX();
						ac[cc].deltaKeyframes[ff].second[1] = qt.getY();
						ac[cc].deltaKeyframes[ff].second[2] = qt.getZ();
						ac[cc].deltaKeyframes[ff].second[3] = qt.getW();
					}
				}

				else if (mac[cc].target_path == "scale")
				{
					ac[cc].typeDelta = 2;

					for (int32 ff = 0; ff < maso.count; ff++)
					{
						float* sca = (float*)&mb[mbv[macso.bufferView].buffer].data.at(offset + ff * stride);
						ac[cc].deltaKeyframes[ff].first = ff;
						ac[cc].deltaKeyframes[ff].second[0] = sca[0];
						ac[cc].deltaKeyframes[ff].second[1] = sca[1];
						ac[cc].deltaKeyframes[ff].second[2] = sca[2];
					}
				}
			}

            Array<float> shapeAnimeWeightArray;									
			for (int32 i= 0;i<ac.size(); i++)									
			{
				Sampler& sa = as[ ac[i].idxSampler ];
				std::pair<int32, float> f0, f1;

				for (int32 kf = 1; kf < sa.interpolateKeyframes.size() ; kf++)   
				{
					f0 = sa.interpolateKeyframes[kf-1];
					f1 = sa.interpolateKeyframes[kf];
					if (f0.second <= frametime && frametime < f1.second ) break;
				}

				float &lowtime = f0.second;
				float &upptime = f1.second;
				const int32 &lowframe = f0.first;
				const int32 &uppframe = f1.first;

				const double mix = (frametime - lowtime) / (upptime - lowtime);

				auto& interpol = mas[ ac[i].idxSampler ].interpolation;
				if      (interpol == "STEP")        gltfInterpolateStep  ( ac[i], lowframe, nodeAniParams );
				else if (interpol == "LINEAR")      gltfInterpolateLinear( ac[i], lowframe, uppframe, mix, nodeAniParams );
				else if (interpol == "CUBICSPLINE") gltfInterpolateSpline( ac[i], lowframe, uppframe, lowtime, upptime, mix, nodeAniParams );

				if (ac[i].idxSampler == -1) continue;

				Array<float>& l = ac[i].deltaKeyframes[lowframe].second;
				Array<float>& u = ac[i].deltaKeyframes[uppframe].second;

				shapeAnimeWeightArray.resize(ac[i].numMorph);
				for (int32 mm = 0; mm < ac[i].numMorph; mm++)
				{
					double weight = l[mm] * (1.0 - mix) + u[mm] * mix;
					shapeAnimeWeightArray[mm] = weight;
				}
			}

			for (int32 nn = 0; nn < gm.nodes.size(); nn++)
			{
				auto& mn = gm.nodes[nn];
				for (int32 cc = 0; cc < mn.children.size(); cc++)
					gltfCalcSkeleton(gm, Mat4x4::Identity(), mn.children[cc], nodeAniParams );
			}

			Array<Array<Mat4x4>> Joints( gm.skins.size() );
			for (int32 nn = 0; nn < gm.nodes.size(); nn++)
			{
				auto& mn = gm.nodes[nn];
				for (int32 cc = 0; cc < mn.children.size(); cc++)
				{
					auto& node = gm.nodes[mn.children[cc]];
					if (node.skin < 0) continue;         

					auto& msns = gm.skins[node.skin];
					auto& ibma = gm.accessors[msns.inverseBindMatrices];
					auto& ibmbv = gm.bufferViews[ibma.bufferView];
					auto  ibmd = gm.buffers[ibmbv.buffer].data.data() + ibma.byteOffset + ibmbv.byteOffset;

					Joints[node.skin].resize( msns.joints.size() );
					for (int32 ii = 0; ii < msns.joints.size(); ii++)
					{
						Mat4x4 ibm = *(Mat4x4*)&ibmd[ii * sizeof(Mat4x4)];
						Mat4x4 matworld = nodeAniParams[ msns.joints[ii] ].matWorld;
						Joints[node.skin][ii] = ibm * matworld;
					}
				}
			}

			for (int32 nn = 0; nn < gm.nodes.size(); nn++)
			{
				for (int32 cc = 0; cc < gm.nodes[nn].children.size(); cc++)
				{
					auto& node = gm.nodes[ gm.nodes[nn].children[cc] ];
					if (node.mesh >= 0)
					{
						uint32 morphidx = 0;                     

						if (cf == 0)
							gltfSetupMorph(node, aniModel.morphMesh); 

						int32 prsize = gltfModel.meshes[node.mesh].primitives.size();
						PrecAnime& precanime = aniModel.precAnimes[animeid];

						for (int32 pp = 0; pp < prsize; pp++)
						{
							auto& pr = gm.meshes[node.mesh].primitives[pp];
							auto& map = gm.accessors[pr.attributes["POSITION"]];

							size_t opos, otex, onor, ojoints, oweights, oidx;
							int32 type_p, type_t, type_n, type_j, type_w, type_i;
							int32 stride_p, stride_t, stride_n, stride_j, stride_w, stride_i;

							auto& bpos = *getBuffer(gm, pr, "POSITION", &opos, &stride_p, &type_p);        
							auto& btex = *getBuffer(gm, pr, "TEXCOORD_0", &otex, &stride_t, &type_t);	  
							auto& bnormal = *getBuffer(gm, pr, "NORMAL", &onor, &stride_n, &type_n);       
							auto& bjoint = *getBuffer(gm, pr, "JOINTS_0", &ojoints, &stride_j, &type_j);   
							auto& bweight = *getBuffer(gm, pr, "WEIGHTS_0", &oweights, &stride_w, &type_w);
							auto& bidx = *getBuffer(gm, pr, &oidx, &stride_i, &type_i);                    

							Array<Vertex3D> vertices(map.count);
							for (int32 vv = 0; vv < map.count; vv++)	
							{
								Vertex3D mv;
								float* pos = (float*)&bpos.data.at(vv * stride_p + opos);
								float* tex = (float*)&btex.data.at(vv * stride_t + otex);
								float* nor = (float*)&bnormal.data.at(vv * stride_n + onor);

								mv.pos = Float3(pos[0], pos[1], pos[2]);
								mv.tex = Float2(tex[0], tex[1]);
								mv.normal = Float3(nor[0], nor[1], nor[2]);

								if (pr.targets.size() && shapeAnimeWeightArray.size())
								{
									for (int32 tt = 0; tt < shapeAnimeWeightArray.size(); tt++)
									{
										if (shapeAnimeWeightArray[tt] == 0) continue;

										size_t opos, onor;
										auto& mtpos = *getBuffer(gm, pr, tt, "POSITION", &opos, &stride_p, &type_p);
										auto& mtnor = *getBuffer(gm, pr, tt, "NORMAL", &onor, &stride_n, &type_n);
										float* spos = (float*)&mtpos.data.at(vv * stride_p + opos);
										float* snor = (float*)&mtnor.data.at(vv * stride_n + onor);
										Float3 shapepos = Float3(spos[0], spos[1], spos[2]);
										Float3 shapenor = Float3(snor[0], snor[1], snor[2]);
										mv.pos += shapepos * shapeAnimeWeightArray[tt];
										mv.normal += shapenor * shapeAnimeWeightArray[tt];
									}
								}

								if (node.skin >= 0)
								{

									uint8* jb = (uint8*)&bjoint.data.at(vv * stride_j + ojoints); 
									uint16* jw = (uint16*)&bjoint.data.at(vv * stride_j + ojoints);
									Word4 j4 = (type_j == TINYGLTF_COMPONENT_TYPE_UNSIGNED_SHORT) ?
														Word4(jw[0], jw[1], jw[2], jw[3]) :
														Word4(jb[0], jb[1], jb[2], jb[3]);

#ifdef USE_VRM
									float* jd = (float*)&bjoint.data.at(vv * stride_j + ojoints);
									Word4 j4 = Word4(jd[0], jd[1], jd[2], jd[3]);
#endif

									float* wf = (float*)&bweight.data.at(vv * stride_w + oweights);
									Float4 w4 = Float4(wf[0], wf[1], wf[2], wf[3]);

									Mat4x4 matskin = w4.x * Joints[node.skin][j4.x] +
										w4.y * Joints[node.skin][j4.y] +
										w4.z * Joints[node.skin][j4.z] +
										w4.w * Joints[node.skin][j4.w];

									if (pr.targets.size() > 0)						
									{
										precanime.Frames[cf].morphMatBuffers.emplace_back(matskin);

										if (aniModel.morphMesh.TexCoord.size() < aniModel.morphMesh.TexCoordCount)
											aniModel.morphMesh.TexCoord.emplace_back(mv.tex);    
									}

									SIMD_Float4 vec4pos = DirectX::XMVector4Transform(SIMD_Float4(mv.pos, 1.0f), matskin);

									Mat4x4 matnor = Mat4x4::Identity();
									if (!(std::abs(matskin.determinant()) < 10e-10))
										matnor = matskin.inverse().transposed();

									mv.pos = vec4pos.xyz() / vec4pos.getW();
									mv.normal = SIMD_Float4{ DirectX::XMVector4Transform(SIMD_Float4(mv.normal, 1.0f), matskin) }.xyz();
								}
								vertices[vv] = mv;							
							}

							MeshData md;
							if (pr.indices > 0)										
							{
								auto& mapi = gm.accessors[pr.indices];
								const uint32 NUMIDX = mapi.count / 3;
								Array<TriangleIndex32> indices(NUMIDX);
								for (int32 ii = 0; ii < NUMIDX; ii += 1)
								{
									TriangleIndex32 idx = TriangleIndex32 ::Zero() ;
									if (mapi.componentType == 5123)
									{
										uint16* ibuf = (uint16*)&bidx.data.at(ii * 3 * 2 + oidx); 
										idx.i0 = ibuf[0]; idx.i1 = ibuf[1]; idx.i2 = ibuf[2];
									}
									else if (mapi.componentType == 5125)
									{
										uint32* ibuf = (uint32*)&bidx.data.at(ii * 3 * 4 + oidx); 
										idx.i0 = ibuf[0]; idx.i1 = ibuf[1]; idx.i2 = ibuf[2];
									}
									indices[ii] = idx;
								}

								md = MeshData{ vertices, indices };
							}

							int32 usetex = 0; 
							Texture tex;
							ColorF col = ColorF(1);

							if (pr.material >= 0)
							{
								auto& nt = gm.materials[pr.material].additionalValues["normalTexture"];  
								int32 idx = -1;
								auto& mmv = gm.materials[pr.material].values;
								auto& bcf = mmv["baseColorFactor"];                                         

								if (mmv.count("baseColorTexture"))
								{
									int32 texidx = mmv["baseColorTexture"].json_double_value["index"];
									idx = gm.textures[texidx].source;
								}

								if (bcf.number_array.size()) col = ColorF(bcf.number_array[0],
																			bcf.number_array[1],
																			bcf.number_array[2],
																			bcf.number_array[3]);
								else                         col = ColorF(1);

								if (idx >= 0 && gm.images.size())
								{
									if (cf == 0)  
									{
										tex = Texture();
										if (gm.images[idx].bufferView >= 0)
										{
											auto& bgfx = gm.bufferViews[gm.images[idx].bufferView];
											auto bimg = &gm.buffers[bgfx.buffer].data.at(bgfx.byteOffset);
											tex = Texture(MemoryReader{ bimg,bgfx.byteLength }, TextureDesc::MippedSRGB);
										}
										else
										{
											auto& mii = gm.images[idx].image;
											tex = Texture(MemoryReader{ (void*)&mii,mii.size() }, TextureDesc::MippedSRGB);
										}
									}
									usetex = 1;
								}
							}

							if (cf == 0)  
							{
								precanime.meshTexs.emplace_back(tex);            
								precanime.meshColors.emplace_back(col);          
							}

							precanime.Frames[cf].MeshDatas.emplace_back(md);			 
							precanime.Frames[cf].Meshes.emplace_back(DynamicMesh{ md }); 
							precanime.Frames[cf].useTex.emplace_back(usetex);			 

						}
					}
				}
			}

			Joints.clear();
			shapeAnimeWeightArray.clear();
			nodeAniParams.clear();

			Float3 vmin = { FLT_MAX,FLT_MAX,FLT_MAX };
			Float3 vmax = { FLT_MIN,FLT_MIN,FLT_MIN };

			for (int32 i = 0; i < aniModel.precAnimes[ animeid ].Frames[cf].MeshDatas.size(); i++)
			{
				for (int32 ii = 0; ii < aniModel.precAnimes[ animeid ].Frames[cf].MeshDatas[i].vertices.size(); ii++)
				{
					Vertex3D& mv = aniModel.precAnimes[ animeid ].Frames[cf].MeshDatas[i].vertices[ii];
					if (vmin.x > mv.pos.x) vmin.x = mv.pos.x;
					if (vmin.y > mv.pos.y) vmin.y = mv.pos.y;
					if (vmin.z > mv.pos.z) vmin.z = mv.pos.z;
					if (vmax.x < mv.pos.x) vmax.x = mv.pos.x;
					if (vmax.y < mv.pos.y) vmax.y = mv.pos.y;
					if (vmax.z < mv.pos.z) vmax.z = mv.pos.z;
				}
			}

			aniModel.precAnimes[ animeid ].Frames[cf].obSize = (vmax - vmin);
			aniModel.precAnimes[ animeid ].Frames[cf].obCenter = vmin + (vmax - vmin)/2;
		}

		for (int32 th = 0; th < tmax; th++)
		{
			auto& as = aniModel.precAnimes[animeid].Samplers[th];
			auto& ac = aniModel.precAnimes[animeid].Channels[th];
			for (int32 ss = 0; ss < mas.size(); ss++) as[ss].interpolateKeyframes.clear();
			for (int32 cc = 0; cc < ac.size(); cc++) ac[cc].deltaKeyframes.clear();
			for (int32 ss = 0; ss < mas.size(); ss++) as[ss].interpolateKeyframes.shrink_to_fit();
			for (int32 cc = 0; cc < ac.size(); cc++) ac[cc].deltaKeyframes.shrink_to_fit();
			ac.clear();
			as.clear();
			ac.shrink_to_fit();
			as.shrink_to_fit();
		}

		for (int32 cf = 0; cf < cycleframe; cf++)
			aniModel.precAnimes[animeid].Frames[cf].MeshDatas.shrink_to_fit();

		aniModel.precAnimes[animeid].Samplers.clear();
		aniModel.precAnimes[animeid].Channels.clear();
		aniModel.precAnimes[animeid].Samplers.shrink_to_fit();
		aniModel.precAnimes[animeid].Channels.shrink_to_fit();

	}

	void gltfInterpolateStep( Channel& ch, int32 lowframe, Array<NodeParam>& _nodeParams )
    {
		Array<float>& v = ch.deltaKeyframes[lowframe].second;
		if		(ch.typeDelta == 1) _nodeParams[ch.idxNode].posePos = Float3{ v[0], v[1], v[2] };
		else if (ch.typeDelta == 3) _nodeParams[ch.idxNode].poseRot = Float4{ v[0], v[1], v[2], v[3] };
		else if	(ch.typeDelta == 2) _nodeParams[ch.idxNode].poseSca = Float3{ v[0], v[1], v[2] };
	}

    void gltfInterpolateLinear( Channel& ch, int32 lowframe, int32 uppframe, float tt, Array<NodeParam>& _nodeParams)
    {
		Array<float>& l = ch.deltaKeyframes[lowframe].second;
		Array<float>& u = ch.deltaKeyframes[uppframe].second;
		if (ch.typeDelta == 1)		
        {
			Float3 low{ l[0], l[1], l[2] };
			Float3 upp{ u[0], u[1], u[2] };
			_nodeParams[ch.idxNode].posePos = low * (1.0 - tt) + upp * tt;
		}

		else if (ch.typeDelta == 3)		
        {
			Float4 low{ l[0], l[1], l[2], l[3] };
			Float4 upp{ u[0], u[1], u[2], l[3] };
            Quaternion lr = Quaternion(low.x, low.y, low.z, low.w);
            Quaternion ur = Quaternion(upp.x, upp.y, upp.z, upp.w);
            Quaternion mx = lr.slerp(ur, tt).normalize();
			_nodeParams[ch.idxNode].poseRot = Float4{ mx.getX(), mx.getY(), mx.getZ(), mx.getW() };
		}

		else if (ch.typeDelta == 2)		
        {
			Float3 low{ l[0], l[1], l[2] };
			Float3 upp{ u[0], u[1], u[2] };
			_nodeParams[ch.idxNode].poseSca = low * (1.0 - tt) + upp * tt;
        }
    }

    template <typename T> T cubicSpline(float tt, T v0, T bb, T v1, T aa)
    {
        const auto t2 = tt * tt;
        const auto t3 = t2 * tt;
        return (2 * t3 - 3 * t2 + 1) * v0 + (t3 - 2 * t2 + tt) * bb + (-2 * t3 + 3 * t2) * v1 + (t3 - t2) * aa;
    }

    void gltfInterpolateSpline( Channel& ch, int32 lowframe,
		                       int32 uppframe, float lowtime, float upptime, float tt, Array<NodeParam>& _nodeParams)
    {
		float delta = upptime - lowtime;

		Array<float>& l = ch.deltaKeyframes[3 * lowframe + 1].second;
		Float4 v0{ l[0],l[1],l[2],l[3] };

		Array<float>& a = ch.deltaKeyframes[3 * uppframe + 0].second;
		Float4 aa = delta * Float4{ a[0],a[1],a[2],a[3] };

		Array<float>& b = ch.deltaKeyframes[3 * lowframe + 2].second;
		Float4 bb = delta * Float4{ b[0],b[1],b[2],b[3] };

		Array<float>& u = ch.deltaKeyframes[3 * uppframe + 1].second;
		Float4 v1{ u[0],u[1],u[2],u[3] };

		if (ch.typeDelta == 1) 
			_nodeParams[ch.idxNode].posePos = cubicSpline(tt, v0, bb, v1, aa).xyz();

		else if (ch.typeDelta == 3) 
        {
            Float4 val = cubicSpline(tt, v0, bb, v1, aa);
            Quaternion qt = Quaternion(val.x, val.y, val.z, val.w).normalize();
			_nodeParams[ch.idxNode].poseRot = qt.toFloat4();
		}

		else if (ch.typeDelta == 2) 
			_nodeParams[ch.idxNode].poseSca = cubicSpline(tt, v0, bb, v1, aa).xyz();
    }

	int32 gltfGetJoint( String jointname )
	{
		for (int32 nn = 0; nn < gltfModel.nodes.size(); nn++)
		{
			auto& msn = gltfModel.nodes[nn];
			if (Unicode::FromUTF8(msn.name) == jointname ) return nn;
		}
		return -1;
	}


	void gltfSetJoint(int32 handle, Float3 trans, Float3 rotate = { 0,0,0 }, Float3 scale = { 1,1,1 } )
	{
		if ( handle < 0 ) return;
		__m128 qrot = XMQuaternionRotationRollPitchYaw(ToRadians(rotate.x), ToRadians(rotate.y), ToRadians(rotate.z)) ;
		Mat4x4 matmodify = Mat4x4::Identity().Scale(scale) * Mat4x4(qrot) *Mat4x4::Identity().Translate(trans);
		nodeParams[handle].matModify = matmodify;
		nodeParams[handle].update = true;
	}

	void gltfSetJoint(int32 handle, Float3 trans, Quaternion rotate = { 0,0,0,1 }, Float3 scale = { 1,1,1 })
	{
		if ( handle < 0 ) return;
		Mat4x4 matmodify = Mat4x4::Identity().Scale(scale) * Mat4x4(rotate) *Mat4x4::Identity().Translate(trans);
		nodeParams[handle].matModify = matmodify;
		nodeParams[handle].update = true;
	}
	Mat4x4 gltfGetMatrix(int32 handle)
	{
		if ( handle < 0 ) return Mat4x4::Identity();
		return nodeParams[handle].matWorld;
	}

    PixieMesh &drawAnime( int32 anime_no = 0,int32 drawframe = NOTUSE, ColorF usrColor=ColorF(0, 0, 0, USE_TEXTURE), int32 istart = NOTUSE, int32 icount = NOTUSE)
    {
		if (Pos.hasNaN() || qRot.hasNaN() || qRot.hasInf()) return *this;

		Rect rectdraw = Rect{ 0,0,camera.getSceneSize() };
        matVP = camera.getViewProj();

		if (!qSpin.isIdentity()) qRot *= qSpin;

        AnimeModel& ani = aniModel;
        __m128 rotate = XMQuaternionRotationRollPitchYaw(ToRadians(eRot.x), ToRadians(eRot.y), ToRadians(eRot.z));
		rotate = qRot * Quaternion(rotate);

		Mat4x4 mrot = Mat4x4(rotate);

        Float3 trans = Pos + rPos;

		Mat = Mat4x4::Identity().Scale(Float3{ -Sca.x,Sca.y,Sca.z }) * mrot * Mat4x4::Identity().Translate(trans);

        PrecAnime& anime = ani.precAnimes[(anime_no == -1) ? 0 : anime_no];
		if (anime.Frames.size() == 0) return *this;

		int32& cf = (drawframe == -1) ? currentFrame : drawframe;

        Frame& frame = anime.Frames[cf];

        uint32 primitiveidx = 0;								
        uint32 tid = 0;

		Array<Vertex3D> morphmv;
		const Array<Array<Vertex3D>>& shapes = ani.morphMesh.ShapeBuffers;

		const uint32 NMORPH = morphTarget.size();
		for (uint32 i = 0; i < frame.Meshes.size(); i++)	
        {
            if ( ani.morphMesh.Targets[i] )
            {
                morphmv = ani.morphMesh.BasisBuffers[primitiveidx]; 

				for (uint32 vv = 0; vv < morphmv.size(); vv++)
                {
                    for (uint32 mm = 0; mm < morphTarget.size(); mm++)
                    {
						float& weight = morphTarget[mm].Weight;
						if (weight == 0.0) continue;
						morphmv[vv].pos += shapes[primitiveidx * NMORPH + mm][vv].pos * weight;
						morphmv[vv].normal += shapes[primitiveidx * NMORPH + mm][vv].normal * weight;
					}

					Mat4x4& matskin = frame.morphMatBuffers[vv];
					SIMD_Float4 vec4pos = DirectX::XMVector4Transform(SIMD_Float4(morphmv[vv].pos, 1.0f), matskin);
					Mat4x4 matnor = matskin.inverse().transposed();

					morphmv[vv].pos = vec4pos.xyz() / vec4pos.getW();
					morphmv[vv].normal = SIMD_Float4{ DirectX::XMVector4Transform(SIMD_Float4(morphmv[vv].normal, 1.0f), matnor) }.xyz();
				}
				frame.Meshes[i].fill(morphmv);				

				primitiveidx++;
			}

			if (istart < 0 )	
			{
				if (usrColor.a >= USE_USERCOLOR)				
					frame.Meshes[i].draw(Mat, ColorF(usrColor.rgb(), 1));

				else			
				{
					if ( usrColor.a == USE_OFFSET_METARIAL )	
						frame.Meshes[i].draw(Mat, anime.meshColors[i] + ColorF(usrColor.rgb(),1) );
					else
						frame.Meshes[i].draw(Mat, anime.meshTexs[i], anime.meshColors[i]);
				}
			}
			else				
			{
				if (usrColor.a >= USE_USERCOLOR)				
					frame.Meshes[i].drawSubset(istart, icount, Mat, anime.meshTexs[tid++]);
				else
				{
					if ( usrColor.a == USE_OFFSET_METARIAL )	
						frame.Meshes[i].drawSubset(istart, icount, Mat, anime.meshColors[i] + ColorF(usrColor.rgb(),1) );
					else
						frame.Meshes[i].drawSubset(istart, icount, Mat, anime.meshTexs[i], anime.meshColors[i]);
				}
            }
        }
		return *this;
    }

	PixieMesh& nextFrame( uint32 anime_no )
    {
        PrecAnime &anime = aniModel.precAnimes[anime_no];
		currentFrame++;
		if ( currentFrame >= (anime.Frames.size()-1) ) currentFrame = 0;
		return *this;
	}




	PixieMesh& drawString( Mat4x4 mat, String text, float kerning = 10,    
		                   ColorF color = Palette::Gray,
		                   int32 istart = 0, float icount = 0)	
	{
		if (0 == noaModel.Meshes.size()) return *this;
		if (mat.value.r->m128_f32[0]==-Math::NaN || Pos.hasNaN() || qRot.hasNaN() || qRot.hasInf()) assert(1);

		const uint8 CODEMAP[256] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13, 0, 0, 
									 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
									 0,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,   
									29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,  
									45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,  
									61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,  
									77,78,79,80,81,82,83,84,85,86,87,88,89,90,91, 
									92,93,94,95,96,97,98,99,100,101,102,103,104,105,  
									106,107, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
									 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
									 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
									 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
									 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
									 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
									 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
									 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

		text += U" ";   
		int32 maxcount = text.size();
		bool isall = false;
		if (icount < 0)  
		{
			icount = abs(icount);
			isall = true;
		}
		else if (icount == 0)
		{
			icount = maxcount;
		}

		icount = (icount > maxcount) ? maxcount : icount;

		int32 start = (istart > maxcount) ? maxcount : istart;
		size_t last = icount;
		float  stroke = icount - (uint32)icount;
		last = ((start + last) > maxcount) ? maxcount : start + last;

		__m128 rot = XMQuaternionRotationRollPitchYaw(ToRadians(eRot.x), ToRadians(eRot.y), ToRadians(eRot.z));

		Mat4x4 mrot;
		if (!qSpin.isIdentity()) mrot = Mat4x4(Quaternion(rot) * qRot * qSpin);
		else mrot = Mat4x4(Quaternion(rot) * qRot);

		mat = mat * Mat4x4::Identity().Scale(Float3{ -Sca.x,Sca.y,Sca.z }) * mrot;

		Float3 pos = Pos;
		Float3 f3 = Mat4x4(qRot).transformPoint(Float3(kerning, 0, 0));
		Float3 lf = Mat4x4(qRot).transformPoint(Float3(0, kerning, 0));
		Float3 cr = pos;

		for (int32 i = start; i <= last; i++)
		{
			auto ascii = text.substr(i, 1)[0];
			if (ascii >= sizeof(CODEMAP)) continue;
			const uint8& code = CODEMAP[ascii];

			if (ascii == ' ')
			{
				pos += f3;
				continue;
			}
			else if (ascii == '\n')
			{
				pos = cr + lf;
				cr = pos;
				continue;
			}

			if (displacementFunc != nullptr)
			{
				Array<Vertex3D>	vertices = noaModel.MeshDatas[code].vertices;	
				Array<TriangleIndex32> indices = noaModel.MeshDatas[code].indices;
				(*displacementFunc)(vertices, indices);
				noaModel.Meshes[code].fill(vertices);							
				noaModel.Meshes[code].fill(indices);
			}

			if (isall)
			{
				auto count = noaModel.Meshes[code].num_triangles();
				count = count * stroke;
				noaModel.Meshes[code].drawSubset(0, count, mat.translated(pos), color);
				pos += f3;
			}
			else
			{
				if (i < (last - 1))
				{
					noaModel.Meshes[code].draw(mat.translated(pos), color);
					pos += f3;
				}
				else noaModel.Meshes[code].drawSubset(0, noaModel.Meshes[code].num_triangles() * stroke, mat.translated(pos), color);
			}
		}
	}





	PixieMesh& drawString(String text, float kerning = 10, float radius = 0, ColorF color = Palette::White,
		int32 istart = 0, float icount = 0 )	
	{
		if (Pos.hasNaN() || qRot.hasNaN() || qRot.hasInf()) return *this;

		if ( 0 == noaModel.Meshes.size()) return *this;
		const uint8 CODEMAP[] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13, 0, 0, 
									 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
									 0,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,   
									29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,  
									45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,  
									61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,  
									77,108,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,  
									93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,  
									 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
									 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
									 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
									 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
									 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
									 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
									 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
									 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

		text += U" ";   
		int32 maxcount = text.size();
		bool isall = false;
		if (icount < 0)  
		{
			icount = abs(icount);
			isall = true;
		}
		else if (icount == 0)
		{
			icount = maxcount;
		}

		icount = (icount > maxcount) ? maxcount : icount;

		int32 start = (istart > maxcount) ? maxcount : istart;
		size_t last = icount;
		float  stroke = icount - (uint32)icount;
		last = ((start + last) > maxcount) ? maxcount : start + last;

		Float3 pos = Pos;
		Mat4x4 mat = Mat4x4::Identity();

		if (radius == 0)	
		{
			__m128 rot = XMQuaternionRotationRollPitchYaw(ToRadians(eRot.x), ToRadians(eRot.y), ToRadians(eRot.z));

			Mat4x4 mrot;
			if (!qSpin.isIdentity()) mrot = Mat4x4(Quaternion(rot) * qRot * qSpin);
			else mrot = Mat4x4(Quaternion(rot) * qRot);

			mat = Mat4x4::Identity().Scale(Float3{ -Sca.x,Sca.y,Sca.z }) * mrot;

			Float3 f3 = Mat4x4(qRot).transformPoint(Float3(kerning, 0, 0));
			Float3 lf = Mat4x4(qRot).transformPoint(Float3(0, kerning, 0));
			Float3 cr = pos;

			for (int32 i = start; i <= last; i++)
			{
				auto ascii = text.substr(i, 1)[0];
				if (ascii >= sizeof(CODEMAP)) continue;
				const uint8& code = CODEMAP[ascii];

				if (ascii == ' ')
				{
					pos += f3;
					continue;
				}
				else if (ascii == '\n')
				{
					pos = cr+lf;
					cr = pos;
					continue;
				}

				if (displacementFunc != nullptr)
				{
					Array<Vertex3D>	vertices = noaModel.MeshDatas[code].vertices;	
					Array<TriangleIndex32> indices = noaModel.MeshDatas[code].indices;
					(*displacementFunc)(vertices, indices);
					noaModel.Meshes[code].fill(vertices);							
					noaModel.Meshes[code].fill(indices);
				}

				if (isall)
				{
					auto count = noaModel.Meshes[code].num_triangles();
					count = count * stroke;
					noaModel.Meshes[code].drawSubset(0, count, mat.translated(pos), color);
					pos += f3;
				}
				else
				{
					if (i < (last - 1))
					{
						noaModel.Meshes[code].draw(mat.translated(pos), color);
						pos += f3;
					}
					else noaModel.Meshes[code].drawSubset(0, noaModel.Meshes[code].num_triangles() * stroke, mat.translated(pos), color);
				}
			}
		}

		else 
		{
			for (uint32 i = start; i <= last; i++)
			{
				auto ascii = text.substr(i, 1)[0];
				if (ascii == ' ' || ascii >= sizeof(CODEMAP)) continue;
				const uint8& code = CODEMAP[ascii];

				Quaternion r = Quaternion::RotateY( ToRadians(-(signed)i * kerning) );							                                      	

				Float3 tt;

				if (!qSpin.isIdentity()) qRot *= qSpin;
				tt = (r * qRot) * Vec3(radius, 0, 0);      　

				mat = Mat4x4::Identity().rotated(r).scaled(Float3{ Sca }).translated(Pos + tt);

				Float3 pp = mat.transformPoint(Pos);
				Float3 up = Float3{0,1,0};     
				Quaternion rotate = camera.getQLookAt(pp, Pos, &up);

				Quaternion er = Quaternion::RollPitchYaw(ToRadians(eRot.x), ToRadians(eRot.y), ToRadians(eRot.z));	
				Quaternion r2 = Quaternion::RollPitchYaw(ToRadians(-90), ToRadians(0), ToRadians(0));				
				mat = mat.Rotate(r2 * rotate).translated(Pos + tt).rotated(er).scaled(Float3{ Sca.x,-Sca.y,Sca.z });	

				if (displacementFunc != nullptr)
				{
					Array<Vertex3D>	vertices = noaModel.MeshDatas[code].vertices;	
					Array<TriangleIndex32> indices = noaModel.MeshDatas[code].indices;
					(*displacementFunc)(vertices, indices);
					noaModel.Meshes[code].fill(vertices);							
					noaModel.Meshes[code].fill(indices);
				}

				if (isall)
				{
					size_t count = noaModel.Meshes[code].num_triangles();
					noaModel.Meshes[code].drawSubset(0, count * stroke, mat, color);
				}
				else
				{
					if (i < (last - 1))
						noaModel.Meshes[code].draw(mat.translated(pos), color);
					else
						noaModel.Meshes[code].drawSubset(0, noaModel.Meshes[code].num_triangles() * stroke, mat, color);
				}
			}
		}
		return *this;
	}
};

