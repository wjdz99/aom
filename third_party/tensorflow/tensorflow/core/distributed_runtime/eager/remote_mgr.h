/* Copyright 2019 The TensorFlow Authors. All Rights Reserved.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
==============================================================================*/

#ifndef TENSORFLOW_CORE_DISTRIBUTED_RUNTIME_EAGER_REMOTE_MGR_H_
#define TENSORFLOW_CORE_DISTRIBUTED_RUNTIME_EAGER_REMOTE_MGR_H_

#include <unordered_map>

#include "tensorflow/core/common_runtime/eager/eager_executor.h"
#include "tensorflow/core/common_runtime/eager/tensor_handle.h"
#include "tensorflow/core/distributed_runtime/eager/remote_tensor_handle.h"
#include "tensorflow/core/platform/mutex.h"

namespace tensorflow {
namespace eager {

// This class manages the states required to setup an eager cluster.
// TODO(fishx): Move remote state from context to this class.
class RemoteMgr {
 public:
  RemoteMgr(bool is_master, EagerContext* ctx)
      : is_master_(is_master), parent_(ctx) {}

  ~RemoteMgr() {
    for (const auto& entry : remote_tensor_handle_map_) {
      entry.second->Unref();
    }
  }

  bool IsMaster() { return is_master_; }

  void AddOperationOutputs(
      const gtl::ArraySlice<tensorflow::TensorHandle*> handles,
      int64 operation_id);

  Status GetTensorHandle(const RemoteTensorHandleInternal& remote_handle,
                         tensorflow::TensorHandle** handle);

  Status DeleteTensorHandle(const RemoteTensorHandleInternal& remote_handle);

  // Helper function to create monotonically increasing ids unique to this
  // context.
  uint64 NextOpId() {
    DCHECK(is_master_);
    mutex_lock l(next_id_mutex_);
    return next_op_id_++;
  }

  // Serialize a TensorHandle(local/remote) to a RemoteTensorHandle.
  Status SerializeRemoteTensorHandle(
      TensorHandle* in, RemoteTensorHandle* out, Device* device,
      const string& device_name,
      const bool serialize_resource_dtype_and_shape = false);

  // Deserialize a RemoteTensorHandle to a TensorHandle(local/remote).
  // The output holds a reference to the TensorHandle.
  Status DeserializeRemoteTensorHandle(const RemoteTensorHandle& in,
                                       TensorHandle** out);

  EagerExecutor& GetOrCreateExecutorForStream(uint64 stream_id);

  void DeleteExecutorForStream(uint64 stream_id);

 protected:
  mutex next_id_mutex_;
  uint64 next_op_id_ GUARDED_BY(next_id_mutex_) = 1;

 private:
  // Returns the op_id and output_num if the given local TensorHandle exists in
  // remote_tensor_handle_map_.
  Status GetRemoteTensorHandle(const tensorflow::TensorHandle* handle,
                               int64* op_id, int32* output_num)
      SHARED_LOCKS_REQUIRED(remote_tensor_handle_mu_);

  Status GetTensorHandleImpl(const RemoteTensorHandleInternal& remote_handle,
                             tensorflow::TensorHandle** handle)
      SHARED_LOCKS_REQUIRED(remote_tensor_handle_mu_);

  bool is_master_;

  using RemoteTensorHandleMap =
      gtl::FlatMap<RemoteTensorHandleInternal, tensorflow::TensorHandle*,
                   RemoteTensorHandleInternalHash,
                   RemoteTensorHandleInternalEquals>;
  mutex remote_tensor_handle_mu_;
  // This map maintains the TensorHandles that are required by remote workers
  // in the cluster. Each map key is generated by the master, so it should be
  // globally unique. This map owns references on the handles it contains.
  RemoteTensorHandleMap remote_tensor_handle_map_
      GUARDED_BY(remote_tensor_handle_mu_);

  EagerContext* parent_;  // not owned.

  mutex executor_map_mu_;
  std::unordered_map<uint64, EagerExecutor> executor_map_
      GUARDED_BY(executor_map_mu_);
};

}  // namespace eager
}  // namespace tensorflow

#endif  // TENSORFLOW_CORE_DISTRIBUTED_RUNTIME_EAGER_REMOTE_MGR_H_
