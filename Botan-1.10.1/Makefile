# Compiler Options
CXX           = g++
LIB_OPT       = -O3 -finline-functions
CHECK_OPT     = -O2
MACH_OPT      = -march=i686 -momit-leaf-frame-pointer
LANG_FLAGS    = -D_REENTRANT -Wno-long-long
WARN_FLAGS    = -W -Wall
SO_OBJ_FLAGS  = -fPIC -fvisibility=hidden
SO_LINK_CMD   = $(CXX) -shared -fPIC -Wl,-soname,$(SONAME)
LINK_TO       = -ldl -lpthread -lrt

# Version Numbers
VERSION       = 1.10.1
SERIES        = 1.10

# Installation Settings
DESTDIR       = /usr/local

BINDIR        = $(DESTDIR)/bin
LIBDIR        = $(DESTDIR)/lib
HEADERDIR     = $(DESTDIR)/include/botan-$(SERIES)/botan
DOCDIR        = $(DESTDIR)/share/doc/botan-$(VERSION)
PKGCONF_DIR   = $(LIBDIR)/pkgconfig

CONFIG_SCRIPT = build/botan-config-1.10
PKGCONFIG     = build/botan-1.10.pc

# Aliases for Common Programs
AR               = ar crs
COPY             = cp
COPY_R           = cp -r
CD               = @cd
ECHO             = @echo
INSTALL_CMD_EXEC = install -m 755
INSTALL_CMD_DATA = install -m 644
LN               = ln -fs
MKDIR            = @mkdir
MKDIR_INSTALL    = @umask 022; mkdir -p -m 755
RANLIB           = true
RM               = @rm -f
RM_R             = @rm -rf

# File Lists
CHECK         = check

HEADERS       = ./src/algo_base/algo_base.h \
                ./src/algo_base/buf_comp.h \
                ./src/algo_base/key_spec.h \
                ./src/algo_base/sym_algo.h \
                ./src/algo_base/symkey.h \
                ./src/algo_factory/algo_factory.h \
                ./src/alloc/allocate.h \
                ./src/alloc/secmem.h \
                ./src/asn1/alg_id.h \
                ./src/asn1/asn1_int.h \
                ./src/asn1/asn1_obj.h \
                ./src/asn1/asn1_oid.h \
                ./src/asn1/asn1_str.h \
                ./src/asn1/ber_dec.h \
                ./src/asn1/der_enc.h \
                ./src/asn1/x509_dn.h \
                ./src/benchmark/benchmark.h \
                ./src/block/aes/aes.h \
                ./src/block/block_cipher.h \
                ./src/block/blowfish/blowfish.h \
                ./src/block/cascade/cascade.h \
                ./src/block/cast/cast128.h \
                ./src/block/cast/cast256.h \
                ./src/block/des/des.h \
                ./src/block/des/desx.h \
                ./src/block/gost_28147/gost_28147.h \
                ./src/block/idea/idea.h \
                ./src/block/kasumi/kasumi.h \
                ./src/block/lion/lion.h \
                ./src/block/lubyrack/lubyrack.h \
                ./src/block/mars/mars.h \
                ./src/block/misty1/misty1.h \
                ./src/block/noekeon/noekeon.h \
                ./src/block/noekeon_simd/noekeon_simd.h \
                ./src/block/rc2/rc2.h \
                ./src/block/rc5/rc5.h \
                ./src/block/rc6/rc6.h \
                ./src/block/safer/safer_sk.h \
                ./src/block/seed/seed.h \
                ./src/block/serpent/serpent.h \
                ./src/block/serpent_simd/serp_simd.h \
                ./src/block/serpent_x86_32/serp_x86_32.h \
                ./src/block/skipjack/skipjack.h \
                ./src/block/square/square.h \
                ./src/block/tea/tea.h \
                ./src/block/twofish/twofish.h \
                ./src/block/xtea/xtea.h \
                ./src/block/xtea_simd/xtea_simd.h \
                ./src/cert/certstore/certstor.h \
                ./src/cert/pkcs10/pkcs10.h \
                ./src/cert/x509ca/x509_ca.h \
                ./src/cert/x509cert/x509_ext.h \
                ./src/cert/x509cert/x509_obj.h \
                ./src/cert/x509cert/x509cert.h \
                ./src/cert/x509crl/crl_ent.h \
                ./src/cert/x509crl/x509_crl.h \
                ./src/cert/x509self/x509self.h \
                ./src/cert/x509store/x509stor.h \
                ./src/checksum/adler32/adler32.h \
                ./src/checksum/crc24/crc24.h \
                ./src/checksum/crc32/crc32.h \
                ./src/codec/base64/base64.h \
                ./src/codec/hex/hex.h \
                ./src/codec/openpgp/openpgp.h \
                ./src/codec/pem/pem.h \
                ./src/constructs/aont/package.h \
                ./src/constructs/cryptobox/cryptobox.h \
                ./src/constructs/fpe_fe1/fpe_fe1.h \
                ./src/constructs/rfc3394/rfc3394.h \
                ./src/constructs/tss/tss.h \
                ./src/engine/dyn_engine/dyn_engine.h \
                ./src/engine/engine.h \
                ./src/entropy/entropy_src.h \
                ./src/filters/basefilt.h \
                ./src/filters/buf_filt.h \
                ./src/filters/codec_filt/b64_filt.h \
                ./src/filters/codec_filt/hex_filt.h \
                ./src/filters/data_snk.h \
                ./src/filters/data_src.h \
                ./src/filters/fd_unix/fd_unix.h \
                ./src/filters/filter.h \
                ./src/filters/filters.h \
                ./src/filters/key_filt.h \
                ./src/filters/modes/cbc/cbc.h \
                ./src/filters/modes/cfb/cfb.h \
                ./src/filters/modes/cts/cts.h \
                ./src/filters/modes/eax/eax.h \
                ./src/filters/modes/ecb/ecb.h \
                ./src/filters/modes/mode_pad/mode_pad.h \
                ./src/filters/modes/xts/xts.h \
                ./src/filters/pipe.h \
                ./src/filters/pk_filts/pk_filts.h \
                ./src/filters/secqueue.h \
                ./src/hash/bmw_512/bmw_512.h \
                ./src/hash/comb4p/comb4p.h \
                ./src/hash/gost_3411/gost_3411.h \
                ./src/hash/has160/has160.h \
                ./src/hash/hash.h \
                ./src/hash/keccak/keccak.h \
                ./src/hash/md2/md2.h \
                ./src/hash/md4/md4.h \
                ./src/hash/md4_x86_32/md4_x86_32.h \
                ./src/hash/md5/md5.h \
                ./src/hash/md5_x86_32/md5_x86_32.h \
                ./src/hash/mdx_hash/mdx_hash.h \
                ./src/hash/par_hash/par_hash.h \
                ./src/hash/rmd128/rmd128.h \
                ./src/hash/rmd160/rmd160.h \
                ./src/hash/sha1/sha160.h \
                ./src/hash/sha1_x86_32/sha1_x86_32.h \
                ./src/hash/sha2_32/sha2_32.h \
                ./src/hash/sha2_64/sha2_64.h \
                ./src/hash/skein/skein_512.h \
                ./src/hash/tiger/tiger.h \
                ./src/hash/whirlpool/whrlpool.h \
                ./src/kdf/kdf.h \
                ./src/kdf/kdf1/kdf1.h \
                ./src/kdf/kdf2/kdf2.h \
                ./src/kdf/mgf1/mgf1.h \
                ./src/kdf/prf_ssl3/prf_ssl3.h \
                ./src/kdf/prf_tls/prf_tls.h \
                ./src/kdf/prf_x942/prf_x942.h \
                ./src/libstate/botan.h \
                ./src/libstate/global_state.h \
                ./src/libstate/init.h \
                ./src/libstate/libstate.h \
                ./src/libstate/look_pk.h \
                ./src/libstate/lookup.h \
                ./src/libstate/oid_lookup/oids.h \
                ./src/libstate/scan_name.h \
                ./src/mac/cbc_mac/cbc_mac.h \
                ./src/mac/cmac/cmac.h \
                ./src/mac/hmac/hmac.h \
                ./src/mac/mac.h \
                ./src/mac/ssl3mac/ssl3_mac.h \
                ./src/mac/x919_mac/x919_mac.h \
                ./src/math/bigint/bigint.h \
                ./src/math/bigint/divide.h \
                ./src/math/ec_gfp/curve_gfp.h \
                ./src/math/ec_gfp/point_gfp.h \
                ./src/math/mp/mp_types.h \
                ./src/math/numbertheory/numthry.h \
                ./src/math/numbertheory/pow_mod.h \
                ./src/math/numbertheory/reducer.h \
                ./src/passhash/bcrypt/bcrypt.h \
                ./src/passhash/passhash9/passhash9.h \
                ./src/pbe/get_pbe.h \
                ./src/pbe/pbe.h \
                ./src/pbe/pbes1/pbes1.h \
                ./src/pbe/pbes2/pbes2.h \
                ./src/pbkdf/pbkdf.h \
                ./src/pbkdf/pbkdf1/pbkdf1.h \
                ./src/pbkdf/pbkdf2/pbkdf2.h \
                ./src/pbkdf/pgps2k/pgp_s2k.h \
                ./src/pk_pad/eme.h \
                ./src/pk_pad/eme1/eme1.h \
                ./src/pk_pad/eme_pkcs/eme_pkcs.h \
                ./src/pk_pad/emsa.h \
                ./src/pk_pad/emsa1/emsa1.h \
                ./src/pk_pad/emsa1_bsi/emsa1_bsi.h \
                ./src/pk_pad/emsa2/emsa2.h \
                ./src/pk_pad/emsa3/emsa3.h \
                ./src/pk_pad/emsa4/emsa4.h \
                ./src/pk_pad/emsa_raw/emsa_raw.h \
                ./src/pk_pad/hash_id/hash_id.h \
                ./src/pubkey/blinding.h \
                ./src/pubkey/dh/dh.h \
                ./src/pubkey/dl_algo/dl_algo.h \
                ./src/pubkey/dl_group/dl_group.h \
                ./src/pubkey/dlies/dlies.h \
                ./src/pubkey/dsa/dsa.h \
                ./src/pubkey/ec_group/ec_group.h \
                ./src/pubkey/ecc_key/ecc_key.h \
                ./src/pubkey/ecdh/ecdh.h \
                ./src/pubkey/ecdsa/ecdsa.h \
                ./src/pubkey/elgamal/elgamal.h \
                ./src/pubkey/gost_3410/gost_3410.h \
                ./src/pubkey/if_algo/if_algo.h \
                ./src/pubkey/keypair/keypair.h \
                ./src/pubkey/nr/nr.h \
                ./src/pubkey/pk_keys.h \
                ./src/pubkey/pk_ops.h \
                ./src/pubkey/pkcs8.h \
                ./src/pubkey/pubkey.h \
                ./src/pubkey/pubkey_enums.h \
                ./src/pubkey/rsa/rsa.h \
                ./src/pubkey/rw/rw.h \
                ./src/pubkey/x509_key.h \
                ./src/rng/auto_rng/auto_rng.h \
                ./src/rng/hmac_rng/hmac_rng.h \
                ./src/rng/randpool/randpool.h \
                ./src/rng/rng.h \
                ./src/rng/x931_rng/x931_rng.h \
                ./src/selftest/selftest.h \
                ./src/ssl/tls_client.h \
                ./src/ssl/tls_connection.h \
                ./src/ssl/tls_exceptn.h \
                ./src/ssl/tls_magic.h \
                ./src/ssl/tls_policy.h \
                ./src/ssl/tls_record.h \
                ./src/ssl/tls_server.h \
                ./src/ssl/tls_session_key.h \
                ./src/ssl/tls_suites.h \
                ./src/stream/arc4/arc4.h \
                ./src/stream/ctr/ctr.h \
                ./src/stream/ofb/ofb.h \
                ./src/stream/salsa20/salsa20.h \
                ./src/stream/stream_cipher.h \
                ./src/stream/turing/turing.h \
                ./src/stream/wid_wake/wid_wake.h \
                ./src/utils/bswap.h \
                ./src/utils/charset.h \
                ./src/utils/cpuid.h \
                ./src/utils/datastor/datastor.h \
                ./src/utils/exceptn.h \
                ./src/utils/get_byte.h \
                ./src/utils/loadstor.h \
                ./src/utils/mem_ops.h \
                ./src/utils/parsing.h \
                ./src/utils/rotate.h \
                ./src/utils/time.h \
                ./src/utils/types.h \
                ./src/utils/ui.h \
                ./src/utils/version.h \
                build/build.h

LIBOBJS       = build/lib/algo_base_symkey.o \
                build/lib/algo_factory.o \
                build/lib/algo_factory_prov_weight.o \
                build/lib/alloc_alloc_mmap_mmap_mem.o \
                build/lib/alloc_mem_pool.o \
                build/lib/alloc_system_alloc_defalloc.o \
                build/lib/asn1_alg_id.o \
                build/lib/asn1_asn1_alt.o \
                build/lib/asn1_asn1_att.o \
                build/lib/asn1_asn1_int.o \
                build/lib/asn1_asn1_oid.o \
                build/lib/asn1_asn1_str.o \
                build/lib/asn1_asn1_tm.o \
                build/lib/asn1_ber_dec.o \
                build/lib/asn1_der_enc.o \
                build/lib/asn1_x509_dn.o \
                build/lib/benchmark.o \
                build/lib/block_aes.o \
                build/lib/block_blowfish_blfs_tab.o \
                build/lib/block_blowfish.o \
                build/lib/block_cascade.o \
                build/lib/block_cast_cast128.o \
                build/lib/block_cast_cast256.o \
                build/lib/block_cast_cast_tab.o \
                build/lib/block_des.o \
                build/lib/block_des_des_tab.o \
                build/lib/block_des_desx.o \
                build/lib/block_gost_28147.o \
                build/lib/block_idea.o \
                build/lib/block_kasumi.o \
                build/lib/block_lion.o \
                build/lib/block_lubyrack.o \
                build/lib/block_mars.o \
                build/lib/block_misty1.o \
                build/lib/block_noekeon.o \
                build/lib/block_noekeon_simd.o \
                build/lib/block_rc2.o \
                build/lib/block_rc5.o \
                build/lib/block_rc6.o \
                build/lib/block_safer_safer_sk.o \
                build/lib/block_seed.o \
                build/lib/block_seed_seed_tab.o \
                build/lib/block_serpent.o \
                build/lib/block_serpent_simd_serp_simd.o \
                build/lib/block_serpent_x86_32_serp_x86_32.o \
                build/lib/block_serpent_x86_32_serp_x86_32_imp.o \
                build/lib/block_skipjack.o \
                build/lib/block_square_sqr_tab.o \
                build/lib/block_square.o \
                build/lib/block_tea.o \
                build/lib/block_twofish_two_tab.o \
                build/lib/block_twofish.o \
                build/lib/block_xtea.o \
                build/lib/block_xtea_simd.o \
                build/lib/cert_certstore_certstor.o \
                build/lib/cert_pkcs10.o \
                build/lib/cert_x509ca_x509_ca.o \
                build/lib/cert_x509cert_x509_ext.o \
                build/lib/cert_x509cert_x509_obj.o \
                build/lib/cert_x509cert.o \
                build/lib/cert_x509crl_crl_ent.o \
                build/lib/cert_x509crl_x509_crl.o \
                build/lib/cert_x509self_x509opt.o \
                build/lib/cert_x509self.o \
                build/lib/cert_x509store_x509stor.o \
                build/lib/checksum_adler32.o \
                build/lib/checksum_crc24.o \
                build/lib/checksum_crc32.o \
                build/lib/codec_base64.o \
                build/lib/codec_hex.o \
                build/lib/codec_openpgp.o \
                build/lib/codec_pem.o \
                build/lib/constructs_aont_package.o \
                build/lib/constructs_cryptobox.o \
                build/lib/constructs_fpe_fe1.o \
                build/lib/constructs_rfc3394.o \
                build/lib/constructs_tss.o \
                build/lib/engine_asm_engine.o \
                build/lib/engine_core_engine_core_modes.o \
                build/lib/engine_core_engine_def_pk_ops.o \
                build/lib/engine_core_engine_def_powm.o \
                build/lib/engine_core_engine_lookup_block.o \
                build/lib/engine_core_engine_lookup_hash.o \
                build/lib/engine_core_engine_lookup_mac.o \
                build/lib/engine_core_engine_lookup_pbkdf.o \
                build/lib/engine_core_engine_lookup_stream.o \
                build/lib/engine_dyn_engine.o \
                build/lib/engine.o \
                build/lib/engine_simd_engine.o \
                build/lib/entropy_dev_random.o \
                build/lib/entropy_egd_es_egd.o \
                build/lib/entropy_hres_timer.o \
                build/lib/entropy_proc_walk_es_ftw.o \
                build/lib/entropy_unix_procs_es_unix.o \
                build/lib/entropy_unix_procs_unix_cmd.o \
                build/lib/entropy_unix_procs_unix_src.o \
                build/lib/filters_algo_filt.o \
                build/lib/filters_basefilt.o \
                build/lib/filters_buf_filt.o \
                build/lib/filters_codec_filt_b64_filt.o \
                build/lib/filters_codec_filt_hex_filt.o \
                build/lib/filters_data_snk.o \
                build/lib/filters_data_src.o \
                build/lib/filters_fd_unix.o \
                build/lib/filters_filter.o \
                build/lib/filters_modes_cbc.o \
                build/lib/filters_modes_cfb.o \
                build/lib/filters_modes_cts.o \
                build/lib/filters_modes_eax.o \
                build/lib/filters_modes_eax_eax_dec.o \
                build/lib/filters_modes_ecb.o \
                build/lib/filters_modes_mode_pad.o \
                build/lib/filters_modes_xts.o \
                build/lib/filters_out_buf.o \
                build/lib/filters_pipe.o \
                build/lib/filters_pipe_io.o \
                build/lib/filters_pipe_rw.o \
                build/lib/filters_pk_filts.o \
                build/lib/filters_secqueue.o \
                build/lib/hash_bmw_512.o \
                build/lib/hash_comb4p.o \
                build/lib/hash_gost_3411.o \
                build/lib/hash_has160.o \
                build/lib/hash_keccak.o \
                build/lib/hash_md2.o \
                build/lib/hash_md4.o \
                build/lib/hash_md4_x86_32.o \
                build/lib/hash_md4_x86_32_md4_x86_32_imp.o \
                build/lib/hash_md5.o \
                build/lib/hash_md5_x86_32.o \
                build/lib/hash_md5_x86_32_md5_x86_32_imp.o \
                build/lib/hash_mdx_hash.o \
                build/lib/hash_par_hash.o \
                build/lib/hash_rmd128.o \
                build/lib/hash_rmd160.o \
                build/lib/hash_sha1_sha160.o \
                build/lib/hash_sha1_x86_32.o \
                build/lib/hash_sha1_x86_32_sha1_x86_32_imp.o \
                build/lib/hash_sha2_32.o \
                build/lib/hash_sha2_64.o \
                build/lib/hash_skein_skein_512.o \
                build/lib/hash_tiger_tig_tab.o \
                build/lib/hash_tiger.o \
                build/lib/hash_whirlpool_whrl_tab.o \
                build/lib/hash_whirlpool_whrlpool.o \
                build/lib/kdf.o \
                build/lib/kdf_kdf1.o \
                build/lib/kdf_kdf2.o \
                build/lib/kdf_mgf1.o \
                build/lib/kdf_prf_ssl3.o \
                build/lib/kdf_prf_tls.o \
                build/lib/kdf_prf_x942.o \
                build/lib/libstate_get_enc.o \
                build/lib/libstate_global_rng.o \
                build/lib/libstate_global_state.o \
                build/lib/libstate_init.o \
                build/lib/libstate.o \
                build/lib/libstate_lookup.o \
                build/lib/libstate_oid_lookup_oids.o \
                build/lib/libstate_policy.o \
                build/lib/libstate_scan_name.o \
                build/lib/mac_cbc_mac.o \
                build/lib/mac_cmac.o \
                build/lib/mac_hmac.o \
                build/lib/mac.o \
                build/lib/mac_ssl3mac_ssl3_mac.o \
                build/lib/mac_x919_mac.o \
                build/lib/math_bigint_big_code.o \
                build/lib/math_bigint_big_io.o \
                build/lib/math_bigint_big_ops2.o \
                build/lib/math_bigint_big_ops3.o \
                build/lib/math_bigint_big_rand.o \
                build/lib/math_bigint.o \
                build/lib/math_bigint_divide.o \
                build/lib/math_ec_gfp_point_gfp.o \
                build/lib/math_mp_mp_asm.o \
                build/lib/math_mp_mp_comba.o \
                build/lib/math_mp_mp_karat.o \
                build/lib/math_mp_mp_misc.o \
                build/lib/math_mp_mp_monty.o \
                build/lib/math_mp_mp_mulop.o \
                build/lib/math_mp_mp_shift.o \
                build/lib/math_numbertheory_dsa_gen.o \
                build/lib/math_numbertheory_jacobi.o \
                build/lib/math_numbertheory_make_prm.o \
                build/lib/math_numbertheory_mp_numth.o \
                build/lib/math_numbertheory_numthry.o \
                build/lib/math_numbertheory_pow_mod.o \
                build/lib/math_numbertheory_powm_fw.o \
                build/lib/math_numbertheory_powm_mnt.o \
                build/lib/math_numbertheory_primes.o \
                build/lib/math_numbertheory_reducer.o \
                build/lib/math_numbertheory_ressol.o \
                build/lib/mutex_noop_mutex_mux_noop.o \
                build/lib/mutex_pthreads_mux_pthr.o \
                build/lib/passhash_bcrypt.o \
                build/lib/passhash_passhash9.o \
                build/lib/pbe_get_pbe.o \
                build/lib/pbe_pbes1.o \
                build/lib/pbe_pbes2.o \
                build/lib/pbkdf_pbkdf1.o \
                build/lib/pbkdf_pbkdf2.o \
                build/lib/pbkdf_pgps2k_pgp_s2k.o \
                build/lib/pk_pad_eme.o \
                build/lib/pk_pad_eme1.o \
                build/lib/pk_pad_eme_pkcs.o \
                build/lib/pk_pad_emsa1.o \
                build/lib/pk_pad_emsa1_bsi.o \
                build/lib/pk_pad_emsa2.o \
                build/lib/pk_pad_emsa3.o \
                build/lib/pk_pad_emsa4.o \
                build/lib/pk_pad_emsa_raw.o \
                build/lib/pk_pad_hash_id.o \
                build/lib/pubkey_blinding.o \
                build/lib/pubkey_dh.o \
                build/lib/pubkey_dl_algo.o \
                build/lib/pubkey_dl_group.o \
                build/lib/pubkey_dlies.o \
                build/lib/pubkey_dsa.o \
                build/lib/pubkey_ec_group.o \
                build/lib/pubkey_ecc_key.o \
                build/lib/pubkey_ecdh.o \
                build/lib/pubkey_ecdsa.o \
                build/lib/pubkey_elgamal.o \
                build/lib/pubkey_gost_3410.o \
                build/lib/pubkey_if_algo.o \
                build/lib/pubkey_keypair.o \
                build/lib/pubkey_nr.o \
                build/lib/pubkey_pk_algs.o \
                build/lib/pubkey_pk_keys.o \
                build/lib/pubkey_pkcs8.o \
                build/lib/pubkey.o \
                build/lib/pubkey_pubkey_enums.o \
                build/lib/pubkey_rsa.o \
                build/lib/pubkey_rw.o \
                build/lib/pubkey_workfactor.o \
                build/lib/pubkey_x509_key.o \
                build/lib/rng_hmac_rng.o \
                build/lib/rng_randpool.o \
                build/lib/rng.o \
                build/lib/rng_x931_rng.o \
                build/lib/selftest.o \
                build/lib/ssl_c_kex.o \
                build/lib/ssl_cert_req.o \
                build/lib/ssl_cert_ver.o \
                build/lib/ssl_finished.o \
                build/lib/ssl_hello.o \
                build/lib/ssl_rec_read.o \
                build/lib/ssl_rec_wri.o \
                build/lib/ssl_s_kex.o \
                build/lib/ssl_tls_client.o \
                build/lib/ssl_tls_handshake_hash.o \
                build/lib/ssl_tls_policy.o \
                build/lib/ssl_tls_server.o \
                build/lib/ssl_tls_session_key.o \
                build/lib/ssl_tls_state.o \
                build/lib/ssl_tls_suites.o \
                build/lib/stream_arc4.o \
                build/lib/stream_ctr.o \
                build/lib/stream_ofb.o \
                build/lib/stream_salsa20.o \
                build/lib/stream_stream_cipher.o \
                build/lib/stream_turing_tur_tab.o \
                build/lib/stream_turing.o \
                build/lib/stream_wid_wake.o \
                build/lib/utils_assert.o \
                build/lib/utils_charset.o \
                build/lib/utils_cpuid.o \
                build/lib/utils_datastor.o \
                build/lib/utils_dyn_load.o \
                build/lib/utils_mlock.o \
                build/lib/utils_parsing.o \
                build/lib/utils_time.o \
                build/lib/utils_ui.o \
                build/lib/utils_version.o

CHECKOBJS     = build/checks/bench.o \
                build/checks/bigint.o \
                build/checks/check.o \
                build/checks/cvc_tests.o \
                build/checks/dolook.o \
                build/checks/ec_tests.o \
                build/checks/ecdh.o \
                build/checks/ecdsa.o \
                build/checks/misc.o \
                build/checks/pk.o \
                build/checks/pk_bench.o \
                build/checks/timer.o \
                build/checks/validate.o \
                build/checks/x509.o

LIB_FLAGS     = $(LIB_OPT) $(MACH_OPT) $(LANG_FLAGS) $(WARN_FLAGS) $(SO_OBJ_FLAGS)
CHECK_FLAGS   = $(CHECK_OPT) $(LANG_FLAGS) $(WARN_FLAGS)

LIBRARIES     = $(STATIC_LIB) $(SHARED_LIB)

LIBNAME       = libbotan
STATIC_LIB    = $(LIBNAME)-$(SERIES).a

SONAME        = $(LIBNAME)-$(SERIES).so.0
SHARED_LIB    = $(SONAME).1

SYMLINK       = $(LIBNAME)-$(SERIES).so

all: $(LIBRARIES)

# Build Commands
build/lib/algo_base_symkey.o: ./src/algo_base/symkey.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/algo_factory.o: ./src/algo_factory/algo_factory.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/algo_factory_prov_weight.o: ./src/algo_factory/prov_weight.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/alloc_alloc_mmap_mmap_mem.o: ./src/alloc/alloc_mmap/mmap_mem.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/alloc_mem_pool.o: ./src/alloc/mem_pool/mem_pool.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/alloc_system_alloc_defalloc.o: ./src/alloc/system_alloc/defalloc.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/asn1_alg_id.o: ./src/asn1/alg_id.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/asn1_asn1_alt.o: ./src/asn1/asn1_alt.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/asn1_asn1_att.o: ./src/asn1/asn1_att.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/asn1_asn1_int.o: ./src/asn1/asn1_int.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/asn1_asn1_oid.o: ./src/asn1/asn1_oid.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/asn1_asn1_str.o: ./src/asn1/asn1_str.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/asn1_asn1_tm.o: ./src/asn1/asn1_tm.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/asn1_ber_dec.o: ./src/asn1/ber_dec.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/asn1_der_enc.o: ./src/asn1/der_enc.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/asn1_x509_dn.o: ./src/asn1/x509_dn.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/benchmark.o: ./src/benchmark/benchmark.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_aes.o: ./src/block/aes/aes.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_blowfish_blfs_tab.o: ./src/block/blowfish/blfs_tab.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_blowfish.o: ./src/block/blowfish/blowfish.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_cascade.o: ./src/block/cascade/cascade.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_cast_cast128.o: ./src/block/cast/cast128.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_cast_cast256.o: ./src/block/cast/cast256.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_cast_cast_tab.o: ./src/block/cast/cast_tab.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_des.o: ./src/block/des/des.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_des_des_tab.o: ./src/block/des/des_tab.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_des_desx.o: ./src/block/des/desx.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_gost_28147.o: ./src/block/gost_28147/gost_28147.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_idea.o: ./src/block/idea/idea.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_kasumi.o: ./src/block/kasumi/kasumi.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_lion.o: ./src/block/lion/lion.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_lubyrack.o: ./src/block/lubyrack/lubyrack.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_mars.o: ./src/block/mars/mars.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_misty1.o: ./src/block/misty1/misty1.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_noekeon.o: ./src/block/noekeon/noekeon.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_noekeon_simd.o: ./src/block/noekeon_simd/noekeon_simd.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_rc2.o: ./src/block/rc2/rc2.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_rc5.o: ./src/block/rc5/rc5.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_rc6.o: ./src/block/rc6/rc6.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_safer_safer_sk.o: ./src/block/safer/safer_sk.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_seed.o: ./src/block/seed/seed.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_seed_seed_tab.o: ./src/block/seed/seed_tab.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_serpent.o: ./src/block/serpent/serpent.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_serpent_simd_serp_simd.o: ./src/block/serpent_simd/serp_simd.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_serpent_x86_32_serp_x86_32.o: ./src/block/serpent_x86_32/serp_x86_32.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_serpent_x86_32_serp_x86_32_imp.o: ./src/block/serpent_x86_32/serp_x86_32_imp.S
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_skipjack.o: ./src/block/skipjack/skipjack.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_square_sqr_tab.o: ./src/block/square/sqr_tab.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_square.o: ./src/block/square/square.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_tea.o: ./src/block/tea/tea.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_twofish_two_tab.o: ./src/block/twofish/two_tab.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_twofish.o: ./src/block/twofish/twofish.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_xtea.o: ./src/block/xtea/xtea.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/block_xtea_simd.o: ./src/block/xtea_simd/xtea_simd.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/cert_certstore_certstor.o: ./src/cert/certstore/certstor.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/cert_pkcs10.o: ./src/cert/pkcs10/pkcs10.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/cert_x509ca_x509_ca.o: ./src/cert/x509ca/x509_ca.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/cert_x509cert_x509_ext.o: ./src/cert/x509cert/x509_ext.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/cert_x509cert_x509_obj.o: ./src/cert/x509cert/x509_obj.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/cert_x509cert.o: ./src/cert/x509cert/x509cert.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/cert_x509crl_crl_ent.o: ./src/cert/x509crl/crl_ent.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/cert_x509crl_x509_crl.o: ./src/cert/x509crl/x509_crl.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/cert_x509self_x509opt.o: ./src/cert/x509self/x509opt.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/cert_x509self.o: ./src/cert/x509self/x509self.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/cert_x509store_x509stor.o: ./src/cert/x509store/x509stor.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/checksum_adler32.o: ./src/checksum/adler32/adler32.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/checksum_crc24.o: ./src/checksum/crc24/crc24.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/checksum_crc32.o: ./src/checksum/crc32/crc32.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/codec_base64.o: ./src/codec/base64/base64.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/codec_hex.o: ./src/codec/hex/hex.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/codec_openpgp.o: ./src/codec/openpgp/openpgp.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/codec_pem.o: ./src/codec/pem/pem.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/constructs_aont_package.o: ./src/constructs/aont/package.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/constructs_cryptobox.o: ./src/constructs/cryptobox/cryptobox.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/constructs_fpe_fe1.o: ./src/constructs/fpe_fe1/fpe_fe1.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/constructs_rfc3394.o: ./src/constructs/rfc3394/rfc3394.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/constructs_tss.o: ./src/constructs/tss/tss.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/engine_asm_engine.o: ./src/engine/asm_engine/asm_engine.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/engine_core_engine_core_modes.o: ./src/engine/core_engine/core_modes.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/engine_core_engine_def_pk_ops.o: ./src/engine/core_engine/def_pk_ops.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/engine_core_engine_def_powm.o: ./src/engine/core_engine/def_powm.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/engine_core_engine_lookup_block.o: ./src/engine/core_engine/lookup_block.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/engine_core_engine_lookup_hash.o: ./src/engine/core_engine/lookup_hash.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/engine_core_engine_lookup_mac.o: ./src/engine/core_engine/lookup_mac.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/engine_core_engine_lookup_pbkdf.o: ./src/engine/core_engine/lookup_pbkdf.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/engine_core_engine_lookup_stream.o: ./src/engine/core_engine/lookup_stream.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/engine_dyn_engine.o: ./src/engine/dyn_engine/dyn_engine.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/engine.o: ./src/engine/engine.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/engine_simd_engine.o: ./src/engine/simd_engine/simd_engine.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/entropy_dev_random.o: ./src/entropy/dev_random/dev_random.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/entropy_egd_es_egd.o: ./src/entropy/egd/es_egd.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/entropy_hres_timer.o: ./src/entropy/hres_timer/hres_timer.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/entropy_proc_walk_es_ftw.o: ./src/entropy/proc_walk/es_ftw.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/entropy_unix_procs_es_unix.o: ./src/entropy/unix_procs/es_unix.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/entropy_unix_procs_unix_cmd.o: ./src/entropy/unix_procs/unix_cmd.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/entropy_unix_procs_unix_src.o: ./src/entropy/unix_procs/unix_src.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/filters_algo_filt.o: ./src/filters/algo_filt.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/filters_basefilt.o: ./src/filters/basefilt.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/filters_buf_filt.o: ./src/filters/buf_filt.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/filters_codec_filt_b64_filt.o: ./src/filters/codec_filt/b64_filt.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/filters_codec_filt_hex_filt.o: ./src/filters/codec_filt/hex_filt.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/filters_data_snk.o: ./src/filters/data_snk.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/filters_data_src.o: ./src/filters/data_src.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/filters_fd_unix.o: ./src/filters/fd_unix/fd_unix.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/filters_filter.o: ./src/filters/filter.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/filters_modes_cbc.o: ./src/filters/modes/cbc/cbc.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/filters_modes_cfb.o: ./src/filters/modes/cfb/cfb.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/filters_modes_cts.o: ./src/filters/modes/cts/cts.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/filters_modes_eax.o: ./src/filters/modes/eax/eax.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/filters_modes_eax_eax_dec.o: ./src/filters/modes/eax/eax_dec.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/filters_modes_ecb.o: ./src/filters/modes/ecb/ecb.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/filters_modes_mode_pad.o: ./src/filters/modes/mode_pad/mode_pad.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/filters_modes_xts.o: ./src/filters/modes/xts/xts.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/filters_out_buf.o: ./src/filters/out_buf.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/filters_pipe.o: ./src/filters/pipe.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/filters_pipe_io.o: ./src/filters/pipe_io.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/filters_pipe_rw.o: ./src/filters/pipe_rw.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/filters_pk_filts.o: ./src/filters/pk_filts/pk_filts.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/filters_secqueue.o: ./src/filters/secqueue.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_bmw_512.o: ./src/hash/bmw_512/bmw_512.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_comb4p.o: ./src/hash/comb4p/comb4p.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_gost_3411.o: ./src/hash/gost_3411/gost_3411.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_has160.o: ./src/hash/has160/has160.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_keccak.o: ./src/hash/keccak/keccak.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_md2.o: ./src/hash/md2/md2.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_md4.o: ./src/hash/md4/md4.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_md4_x86_32.o: ./src/hash/md4_x86_32/md4_x86_32.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_md4_x86_32_md4_x86_32_imp.o: ./src/hash/md4_x86_32/md4_x86_32_imp.S
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_md5.o: ./src/hash/md5/md5.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_md5_x86_32.o: ./src/hash/md5_x86_32/md5_x86_32.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_md5_x86_32_md5_x86_32_imp.o: ./src/hash/md5_x86_32/md5_x86_32_imp.S
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_mdx_hash.o: ./src/hash/mdx_hash/mdx_hash.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_par_hash.o: ./src/hash/par_hash/par_hash.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_rmd128.o: ./src/hash/rmd128/rmd128.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_rmd160.o: ./src/hash/rmd160/rmd160.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_sha1_sha160.o: ./src/hash/sha1/sha160.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_sha1_x86_32.o: ./src/hash/sha1_x86_32/sha1_x86_32.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_sha1_x86_32_sha1_x86_32_imp.o: ./src/hash/sha1_x86_32/sha1_x86_32_imp.S
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_sha2_32.o: ./src/hash/sha2_32/sha2_32.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_sha2_64.o: ./src/hash/sha2_64/sha2_64.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_skein_skein_512.o: ./src/hash/skein/skein_512.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_tiger_tig_tab.o: ./src/hash/tiger/tig_tab.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_tiger.o: ./src/hash/tiger/tiger.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_whirlpool_whrl_tab.o: ./src/hash/whirlpool/whrl_tab.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/hash_whirlpool_whrlpool.o: ./src/hash/whirlpool/whrlpool.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/kdf.o: ./src/kdf/kdf.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/kdf_kdf1.o: ./src/kdf/kdf1/kdf1.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/kdf_kdf2.o: ./src/kdf/kdf2/kdf2.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/kdf_mgf1.o: ./src/kdf/mgf1/mgf1.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/kdf_prf_ssl3.o: ./src/kdf/prf_ssl3/prf_ssl3.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/kdf_prf_tls.o: ./src/kdf/prf_tls/prf_tls.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/kdf_prf_x942.o: ./src/kdf/prf_x942/prf_x942.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/libstate_get_enc.o: ./src/libstate/get_enc.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/libstate_global_rng.o: ./src/libstate/global_rng.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/libstate_global_state.o: ./src/libstate/global_state.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/libstate_init.o: ./src/libstate/init.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/libstate.o: ./src/libstate/libstate.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/libstate_lookup.o: ./src/libstate/lookup.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/libstate_oid_lookup_oids.o: ./src/libstate/oid_lookup/oids.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/libstate_policy.o: ./src/libstate/policy.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/libstate_scan_name.o: ./src/libstate/scan_name.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/mac_cbc_mac.o: ./src/mac/cbc_mac/cbc_mac.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/mac_cmac.o: ./src/mac/cmac/cmac.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/mac_hmac.o: ./src/mac/hmac/hmac.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/mac.o: ./src/mac/mac.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/mac_ssl3mac_ssl3_mac.o: ./src/mac/ssl3mac/ssl3_mac.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/mac_x919_mac.o: ./src/mac/x919_mac/x919_mac.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_bigint_big_code.o: ./src/math/bigint/big_code.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_bigint_big_io.o: ./src/math/bigint/big_io.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_bigint_big_ops2.o: ./src/math/bigint/big_ops2.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_bigint_big_ops3.o: ./src/math/bigint/big_ops3.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_bigint_big_rand.o: ./src/math/bigint/big_rand.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_bigint.o: ./src/math/bigint/bigint.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_bigint_divide.o: ./src/math/bigint/divide.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_ec_gfp_point_gfp.o: ./src/math/ec_gfp/point_gfp.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_mp_mp_asm.o: ./src/math/mp/mp_asm.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_mp_mp_comba.o: ./src/math/mp/mp_comba.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_mp_mp_karat.o: ./src/math/mp/mp_karat.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_mp_mp_misc.o: ./src/math/mp/mp_misc.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_mp_mp_monty.o: ./src/math/mp/mp_monty.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_mp_mp_mulop.o: ./src/math/mp/mp_mulop.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_mp_mp_shift.o: ./src/math/mp/mp_shift.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_numbertheory_dsa_gen.o: ./src/math/numbertheory/dsa_gen.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_numbertheory_jacobi.o: ./src/math/numbertheory/jacobi.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_numbertheory_make_prm.o: ./src/math/numbertheory/make_prm.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_numbertheory_mp_numth.o: ./src/math/numbertheory/mp_numth.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_numbertheory_numthry.o: ./src/math/numbertheory/numthry.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_numbertheory_pow_mod.o: ./src/math/numbertheory/pow_mod.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_numbertheory_powm_fw.o: ./src/math/numbertheory/powm_fw.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_numbertheory_powm_mnt.o: ./src/math/numbertheory/powm_mnt.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_numbertheory_primes.o: ./src/math/numbertheory/primes.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_numbertheory_reducer.o: ./src/math/numbertheory/reducer.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/math_numbertheory_ressol.o: ./src/math/numbertheory/ressol.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/mutex_noop_mutex_mux_noop.o: ./src/mutex/noop_mutex/mux_noop.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/mutex_pthreads_mux_pthr.o: ./src/mutex/pthreads/mux_pthr.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/passhash_bcrypt.o: ./src/passhash/bcrypt/bcrypt.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/passhash_passhash9.o: ./src/passhash/passhash9/passhash9.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pbe_get_pbe.o: ./src/pbe/get_pbe.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pbe_pbes1.o: ./src/pbe/pbes1/pbes1.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pbe_pbes2.o: ./src/pbe/pbes2/pbes2.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pbkdf_pbkdf1.o: ./src/pbkdf/pbkdf1/pbkdf1.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pbkdf_pbkdf2.o: ./src/pbkdf/pbkdf2/pbkdf2.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pbkdf_pgps2k_pgp_s2k.o: ./src/pbkdf/pgps2k/pgp_s2k.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pk_pad_eme.o: ./src/pk_pad/eme.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pk_pad_eme1.o: ./src/pk_pad/eme1/eme1.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pk_pad_eme_pkcs.o: ./src/pk_pad/eme_pkcs/eme_pkcs.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pk_pad_emsa1.o: ./src/pk_pad/emsa1/emsa1.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pk_pad_emsa1_bsi.o: ./src/pk_pad/emsa1_bsi/emsa1_bsi.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pk_pad_emsa2.o: ./src/pk_pad/emsa2/emsa2.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pk_pad_emsa3.o: ./src/pk_pad/emsa3/emsa3.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pk_pad_emsa4.o: ./src/pk_pad/emsa4/emsa4.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pk_pad_emsa_raw.o: ./src/pk_pad/emsa_raw/emsa_raw.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pk_pad_hash_id.o: ./src/pk_pad/hash_id/hash_id.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey_blinding.o: ./src/pubkey/blinding.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey_dh.o: ./src/pubkey/dh/dh.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey_dl_algo.o: ./src/pubkey/dl_algo/dl_algo.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey_dl_group.o: ./src/pubkey/dl_group/dl_group.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey_dlies.o: ./src/pubkey/dlies/dlies.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey_dsa.o: ./src/pubkey/dsa/dsa.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey_ec_group.o: ./src/pubkey/ec_group/ec_group.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey_ecc_key.o: ./src/pubkey/ecc_key/ecc_key.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey_ecdh.o: ./src/pubkey/ecdh/ecdh.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey_ecdsa.o: ./src/pubkey/ecdsa/ecdsa.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey_elgamal.o: ./src/pubkey/elgamal/elgamal.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey_gost_3410.o: ./src/pubkey/gost_3410/gost_3410.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey_if_algo.o: ./src/pubkey/if_algo/if_algo.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey_keypair.o: ./src/pubkey/keypair/keypair.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey_nr.o: ./src/pubkey/nr/nr.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey_pk_algs.o: ./src/pubkey/pk_algs.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey_pk_keys.o: ./src/pubkey/pk_keys.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey_pkcs8.o: ./src/pubkey/pkcs8.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey.o: ./src/pubkey/pubkey.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey_pubkey_enums.o: ./src/pubkey/pubkey_enums.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey_rsa.o: ./src/pubkey/rsa/rsa.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey_rw.o: ./src/pubkey/rw/rw.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey_workfactor.o: ./src/pubkey/workfactor.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/pubkey_x509_key.o: ./src/pubkey/x509_key.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/rng_hmac_rng.o: ./src/rng/hmac_rng/hmac_rng.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/rng_randpool.o: ./src/rng/randpool/randpool.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/rng.o: ./src/rng/rng.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/rng_x931_rng.o: ./src/rng/x931_rng/x931_rng.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/selftest.o: ./src/selftest/selftest.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/ssl_c_kex.o: ./src/ssl/c_kex.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/ssl_cert_req.o: ./src/ssl/cert_req.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/ssl_cert_ver.o: ./src/ssl/cert_ver.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/ssl_finished.o: ./src/ssl/finished.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/ssl_hello.o: ./src/ssl/hello.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/ssl_rec_read.o: ./src/ssl/rec_read.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/ssl_rec_wri.o: ./src/ssl/rec_wri.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/ssl_s_kex.o: ./src/ssl/s_kex.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/ssl_tls_client.o: ./src/ssl/tls_client.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/ssl_tls_handshake_hash.o: ./src/ssl/tls_handshake_hash.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/ssl_tls_policy.o: ./src/ssl/tls_policy.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/ssl_tls_server.o: ./src/ssl/tls_server.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/ssl_tls_session_key.o: ./src/ssl/tls_session_key.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/ssl_tls_state.o: ./src/ssl/tls_state.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/ssl_tls_suites.o: ./src/ssl/tls_suites.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/stream_arc4.o: ./src/stream/arc4/arc4.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/stream_ctr.o: ./src/stream/ctr/ctr.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/stream_ofb.o: ./src/stream/ofb/ofb.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/stream_salsa20.o: ./src/stream/salsa20/salsa20.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/stream_stream_cipher.o: ./src/stream/stream_cipher.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/stream_turing_tur_tab.o: ./src/stream/turing/tur_tab.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/stream_turing.o: ./src/stream/turing/turing.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/stream_wid_wake.o: ./src/stream/wid_wake/wid_wake.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/utils_assert.o: ./src/utils/assert.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/utils_charset.o: ./src/utils/charset.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/utils_cpuid.o: ./src/utils/cpuid.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/utils_datastor.o: ./src/utils/datastor/datastor.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/utils_dyn_load.o: ./src/utils/dyn_load/dyn_load.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/utils_mlock.o: ./src/utils/mlock.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/utils_parsing.o: ./src/utils/parsing.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/utils_time.o: ./src/utils/time.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/utils_ui.o: ./src/utils/ui.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@

build/lib/utils_version.o: ./src/utils/version.cpp
	$(CXX) -Ibuild/include $(LIB_FLAGS) -c $? -o $@


build/checks/bench.o: ./checks/bench.cpp
	$(CXX) -Ibuild/include $(CHECK_FLAGS) -c $? -o $@

build/checks/bigint.o: ./checks/bigint.cpp
	$(CXX) -Ibuild/include $(CHECK_FLAGS) -c $? -o $@

build/checks/check.o: ./checks/check.cpp
	$(CXX) -Ibuild/include $(CHECK_FLAGS) -c $? -o $@

build/checks/cvc_tests.o: ./checks/cvc_tests.cpp
	$(CXX) -Ibuild/include $(CHECK_FLAGS) -c $? -o $@

build/checks/dolook.o: ./checks/dolook.cpp
	$(CXX) -Ibuild/include $(CHECK_FLAGS) -c $? -o $@

build/checks/ec_tests.o: ./checks/ec_tests.cpp
	$(CXX) -Ibuild/include $(CHECK_FLAGS) -c $? -o $@

build/checks/ecdh.o: ./checks/ecdh.cpp
	$(CXX) -Ibuild/include $(CHECK_FLAGS) -c $? -o $@

build/checks/ecdsa.o: ./checks/ecdsa.cpp
	$(CXX) -Ibuild/include $(CHECK_FLAGS) -c $? -o $@

build/checks/misc.o: ./checks/misc.cpp
	$(CXX) -Ibuild/include $(CHECK_FLAGS) -c $? -o $@

build/checks/pk.o: ./checks/pk.cpp
	$(CXX) -Ibuild/include $(CHECK_FLAGS) -c $? -o $@

build/checks/pk_bench.o: ./checks/pk_bench.cpp
	$(CXX) -Ibuild/include $(CHECK_FLAGS) -c $? -o $@

build/checks/timer.o: ./checks/timer.cpp
	$(CXX) -Ibuild/include $(CHECK_FLAGS) -c $? -o $@

build/checks/validate.o: ./checks/validate.cpp
	$(CXX) -Ibuild/include $(CHECK_FLAGS) -c $? -o $@

build/checks/x509.o: ./checks/x509.cpp
	$(CXX) -Ibuild/include $(CHECK_FLAGS) -c $? -o $@


# Link Commands
$(CHECK): $(LIBRARIES) $(CHECKOBJS)
	$(CXX) $(LDFLAGS) $(CHECKOBJS) $(SHARED_LIB) $(LINK_TO) -o $(CHECK)

$(STATIC_LIB): $(LIBOBJS)
	$(RM) $(STATIC_LIB)
	$(AR) $(STATIC_LIB) $(LIBOBJS)
	$(RANLIB) $(STATIC_LIB)

$(SHARED_LIB): $(LIBOBJS)
	$(SO_LINK_CMD) $(LDFLAGS) $(LIBOBJS) -o $(SHARED_LIB) $(LINK_TO)
	$(LN) $(SHARED_LIB) $(SONAME)
	$(LN) $(SHARED_LIB) $(SYMLINK)

# Fake Targets
.PHONY = docs clean distclean install static shared

static: $(STATIC_LIB)

shared: $(SHARED_LIB)

docs:
	$(COPY) readme.txt build/docs
	$(COPY) doc/*.txt build/docs/manual

clean:
	$(RM_R) build/lib/* build/checks/*
	$(RM) $(LIBRARIES) $(SYMLINK) $(SONAME) $(CHECK)

distclean: clean
	$(RM_R) build
	$(RM) Makefile* $(CONFIG_SCRIPT) $(PKGCONFIG)
	$(RM) botan_all.cpp botan_all.h

install: $(LIBRARIES) docs
	$(ECHO) "Installing Botan into $(DESTDIR)... "
	$(MKDIR_INSTALL) $(DOCDIR)
	$(COPY_R) build/docs/* $(DOCDIR)

	$(MKDIR_INSTALL) $(HEADERDIR)
	for i in $(HEADERS); do \
	   $(INSTALL_CMD_DATA) $$i $(HEADERDIR); \
	 done

	$(MKDIR_INSTALL) $(LIBDIR)
	$(INSTALL_CMD_DATA) $(STATIC_LIB) $(LIBDIR)
	$(INSTALL_CMD_EXEC) $(SHARED_LIB) $(LIBDIR)
	$(CD) $(LIBDIR); $(LN) $(SHARED_LIB) $(SYMLINK)
	$(CD) $(LIBDIR); $(LN) $(SHARED_LIB) $(SONAME)

	$(MKDIR_INSTALL) $(BINDIR)
	$(INSTALL_CMD_EXEC) $(CONFIG_SCRIPT) $(BINDIR)

	$(MKDIR_INSTALL) $(PKGCONF_DIR)
	$(INSTALL_CMD_DATA) $(PKGCONFIG) $(PKGCONF_DIR)
